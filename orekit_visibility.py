from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import csv
import math


@dataclass(frozen=True)
class GroundSensor:
    sensor_id: str
    latitude_deg: float
    longitude_deg: float
    altitude_m: float
    max_solar_elevation_deg: float


def _read_3le(path: str | Path) -> List[Tuple[str, str, str]]:
    lines = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) % 3 != 0:
        raise ValueError(f"Le fichier 3le {path} doit contenir des triplets nom/ligne1/ligne2.")

    triplets: List[Tuple[str, str, str]] = []
    for i in range(0, len(lines), 3):
        name, l1, l2 = lines[i], lines[i + 1], lines[i + 2]
        if not l1.startswith("1 ") or not l2.startswith("2 "):
            raise ValueError(f"Triplet invalide dans {path} pour l'objet '{name}' (format TLE attendu).")
        triplets.append((name, l1, l2))
    return triplets


def _read_ground_sensors_csv(path: str | Path) -> List[GroundSensor]:
    sensors: List[GroundSensor] = []
    with Path(path).open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"sensor_id", "latitude_deg", "longitude_deg", "altitude_m", "max_solar_elevation_deg"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(
                "Le CSV capteurs sol doit contenir les colonnes: "
                "sensor_id, latitude_deg, longitude_deg, altitude_m, max_solar_elevation_deg"
            )
        for row in reader:
            sensors.append(
                GroundSensor(
                    sensor_id=row["sensor_id"],
                    latitude_deg=float(row["latitude_deg"]),
                    longitude_deg=float(row["longitude_deg"]),
                    altitude_m=float(row["altitude_m"]),
                    max_solar_elevation_deg=float(row["max_solar_elevation_deg"]),
                )
            )
    return sensors


def _read_masks_csv(path: str | Path) -> Dict[str, List[Tuple[float, float]]]:
    masks: Dict[str, List[Tuple[float, float]]] = {}
    with Path(path).open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required = {"sensor_id", "azimuth_deg", "min_elevation_deg"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError("Le CSV masques doit contenir les colonnes: sensor_id, azimuth_deg, min_elevation_deg")
        for row in reader:
            sid = row["sensor_id"]
            masks.setdefault(sid, []).append((float(row["azimuth_deg"]) % 360.0, float(row["min_elevation_deg"])))

    for sid, points in masks.items():
        points.sort(key=lambda x: x[0])
        if points and points[0][0] != points[-1][0]:
            points.append((points[0][0] + 360.0, points[0][1]))
        masks[sid] = points
    return masks


def _interp_mask_min_elevation(mask_points: List[Tuple[float, float]], azimuth_deg: float) -> float:
    if not mask_points:
        return -90.0
    az = azimuth_deg % 360.0
    if az < mask_points[0][0]:
        az += 360.0
    for i in range(len(mask_points) - 1):
        a0, e0 = mask_points[i]
        a1, e1 = mask_points[i + 1]
        if a0 <= az <= a1:
            if a1 == a0:
                return e0
            t = (az - a0) / (a1 - a0)
            return e0 + t * (e1 - e0)
    return mask_points[-1][1]


def compute_visibility_metrics(
    target_3le_path: str | Path,
    instants_utc: Sequence[str],
    propagator_type: str,
    output_csv_path: str | Path,
    ground_sensors_csv_path: Optional[str | Path] = None,
    ground_masks_csv_path: Optional[str | Path] = None,
    space_sensors_3le_path: Optional[str | Path] = None,
    orekit_data_path: Optional[str | Path] = None,
) -> Path:
    """
    Calcule les visibilités via Orekit et écrit un CSV de métriques.

    Sortie: sensor_id, object_id, instant_utc, apparent_rate_deg_s, distance_km, phase_angle_deg, moon_angle_deg
    """
    try:
        import orekit
        from orekit.pyhelpers import setup_orekit_curdir, setup_orekit_data
        from org.hipparchus.geometry.euclidean.threed import Vector3D
        from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
        from org.orekit.bodies import CelestialBodyFactory, GeodeticPoint, OneAxisEllipsoid
        from org.orekit.frames import FramesFactory, TopocentricFrame
        from org.orekit.orbits import KeplerianOrbit, PositionAngleType
        from org.orekit.propagation.analytical import KeplerianPropagator
        from org.orekit.propagation.analytical.tle import SGP4, TLE
        from org.orekit.propagation.numerical import NumericalPropagator
        from org.orekit.time import AbsoluteDate, TimeScalesFactory
        from org.orekit.utils import Constants, IERSConventions, PVCoordinates
    except Exception as exc:
        raise RuntimeError("Orekit (wrapper Python) n'est pas disponible dans l'environnement.") from exc

    orekit.initVM()
    if orekit_data_path:
        setup_orekit_data(str(orekit_data_path))
    else:
        setup_orekit_curdir()

    targets = _read_3le(target_3le_path)
    ground_sensors = _read_ground_sensors_csv(ground_sensors_csv_path) if ground_sensors_csv_path else []
    masks = _read_masks_csv(ground_masks_csv_path) if ground_masks_csv_path else {}
    space_sensors = _read_3le(space_sensors_3le_path) if space_sensors_3le_path else []

    utc = TimeScalesFactory.getUTC()
    eme2000 = FramesFactory.getEME2000()
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING, itrf)
    sun = CelestialBodyFactory.getSun()
    moon = CelestialBodyFactory.getMoon()

    def make_tle_propagator(l1: str, l2: str):
        return SGP4.selectExtrapolator(TLE(l1, l2))

    ptype = propagator_type.lower().strip()

    def make_propagator_from_tle(l1: str, l2: str):
        if ptype == "sgp4":
            return make_tle_propagator(l1, l2)
        if ptype == "keplerien":
            tle_prop = make_tle_propagator(l1, l2)
            return KeplerianPropagator(KeplerianOrbit(tle_prop.getInitialState().getOrbit()))
        if ptype in {"numérique", "numerique"}:
            tle_prop = make_tle_propagator(l1, l2)
            state = tle_prop.getInitialState()
            integ = DormandPrince853Integrator(0.1, 120.0, 1.0e-9, 1.0e-9)
            num = NumericalPropagator(integ)
            num.setOrbitType(state.getOrbit().getType())
            num.setPositionAngleType(PositionAngleType.MEAN)
            num.setInitialState(state)
            return num
        raise ValueError("propagator_type doit être 'keplerien', 'sgp4' ou 'numérique'.")

    target_props = {name: make_propagator_from_tle(l1, l2) for name, l1, l2 in targets}
    space_props = {name: make_propagator_from_tle(l1, l2) for name, l1, l2 in space_sensors}

    ground_frames: Dict[str, TopocentricFrame] = {}
    for gs in ground_sensors:
        gp = GeodeticPoint(math.radians(gs.latitude_deg), math.radians(gs.longitude_deg), gs.altitude_m)
        ground_frames[gs.sensor_id] = TopocentricFrame(earth, gp, gs.sensor_id)

    out_path = Path(output_csv_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow([
            "sensor_id",
            "object_id",
            "instant_utc",
            "apparent_rate_deg_s",
            "distance_km",
            "phase_angle_deg",
            "moon_angle_deg",
        ])

        for instant in instants_utc:
            date = AbsoluteDate(instant, utc)
            sun_pv = sun.getPVCoordinates(date, eme2000)
            moon_pv = moon.getPVCoordinates(date, eme2000)

            target_pv: Dict[str, PVCoordinates] = {
                oid: prop.propagate(date).getPVCoordinates(eme2000) for oid, prop in target_props.items()
            }

            for gs in ground_sensors:
                frame = ground_frames[gs.sensor_id]
                solar_el = math.degrees(frame.getElevation(sun_pv.getPosition(), eme2000, date))
                if solar_el > gs.max_solar_elevation_deg:
                    continue

                sensor_pv = frame.getPVCoordinates(date, eme2000)
                sensor_pos = sensor_pv.getPosition()
                moon_dir_sensor = moon_pv.getPosition().subtract(sensor_pos)
                mask = masks.get(gs.sensor_id, [])

                for oid, obj_pv in target_pv.items():
                    az = math.degrees(frame.getAzimuth(obj_pv.getPosition(), eme2000, date)) % 360.0
                    el = math.degrees(frame.getElevation(obj_pv.getPosition(), eme2000, date))
                    min_el = _interp_mask_min_elevation(mask, az) if mask else 0.0
                    if el < min_el:
                        continue

                    rel = obj_pv.getPosition().subtract(sensor_pos)
                    rel_v = obj_pv.getVelocity().subtract(sensor_pv.getVelocity())
                    dist_km = rel.getNorm() / 1000.0
                    rate = math.degrees(Vector3D.crossProduct(rel, rel_v).getNorm() / rel.getNormSq())

                    v_sun = sun_pv.getPosition().subtract(obj_pv.getPosition())
                    v_sensor = sensor_pos.subtract(obj_pv.getPosition())
                    phase = math.degrees(Vector3D.angle(v_sun, v_sensor))
                    moon_angle = math.degrees(Vector3D.angle(rel, moon_dir_sensor))

                    writer.writerow([gs.sensor_id, oid, instant, rate, dist_km, phase, moon_angle])

            for sid, sprop in space_props.items():
                sensor_pv = sprop.propagate(date).getPVCoordinates(eme2000)
                sensor_pos = sensor_pv.getPosition()
                moon_dir_sensor = moon_pv.getPosition().subtract(sensor_pos)
                earth_center_sensor = sensor_pos.negate()
                earth_ang_rad = math.asin(Constants.WGS84_EARTH_EQUATORIAL_RADIUS / sensor_pos.getNorm())

                for oid, obj_pv in target_pv.items():
                    rel = obj_pv.getPosition().subtract(sensor_pos)
                    los_to_earth = Vector3D.angle(rel, earth_center_sensor)
                    if los_to_earth < earth_ang_rad:
                        continue

                    rel_v = obj_pv.getVelocity().subtract(sensor_pv.getVelocity())
                    dist_km = rel.getNorm() / 1000.0
                    rate = math.degrees(Vector3D.crossProduct(rel, rel_v).getNorm() / rel.getNormSq())
                    v_sun = sun_pv.getPosition().subtract(obj_pv.getPosition())
                    v_sensor = sensor_pos.subtract(obj_pv.getPosition())
                    phase = math.degrees(Vector3D.angle(v_sun, v_sensor))
                    moon_angle = math.degrees(Vector3D.angle(rel, moon_dir_sensor))

                    writer.writerow([sid, oid, instant, rate, dist_km, phase, moon_angle])

    return out_path
