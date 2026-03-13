import tempfile
import unittest
from pathlib import Path

from orekit_visibility import (
    _interp_mask_min_elevation,
    _read_3le,
    _read_ground_sensors_csv,
    _read_masks_csv,
    compute_visibility_metrics,
)


class TestOrekitVisibilityHelpers(unittest.TestCase):
    def test_parse_3le_ok(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "targets.3le"
            p.write_text(
                "OBJ\n"
                "1 00005U 58002B   20062.51782528  .00000023  00000-0  28098-4 0  9990\n"
                "2 00005  34.2682 331.5174 1849677 331.7664  19.3264 10.82419157413667\n",
                encoding="utf-8",
            )
            data = _read_3le(p)
            self.assertEqual(len(data), 1)
            self.assertEqual(data[0][0], "OBJ")

    def test_parse_3le_bad_format(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "targets.3le"
            p.write_text("OBJ\nNOT_LINE1\nNOT_LINE2\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                _read_3le(p)

    def test_mask_interp(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "mask.csv"
            p.write_text(
                "sensor_id,azimuth_deg,min_elevation_deg\n"
                "s1,0,10\n"
                "s1,180,20\n",
                encoding="utf-8",
            )
            masks = _read_masks_csv(p)
            self.assertAlmostEqual(_interp_mask_min_elevation(masks["s1"], 90), 15.0)

    def test_ground_sensor_csv(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "ground.csv"
            p.write_text(
                "sensor_id,latitude_deg,longitude_deg,altitude_m,max_solar_elevation_deg\n"
                "rennes,48.1173,-1.6778,45,-6\n",
                encoding="utf-8",
            )
            sensors = _read_ground_sensors_csv(p)
            self.assertEqual(sensors[0].sensor_id, "rennes")


class TestSimulatedRun(unittest.TestCase):
    def test_compute_visibility_simulated_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "targets.3le").write_text(
                "ISS (ZARYA)\n"
                "1 25544U 98067A   25067.51041667  .00016717  00000+0  10270-3 0  9997\n"
                "2 25544  51.6416  35.1023 0004513 165.9410 315.5184 15.50012345678901\n",
                encoding="utf-8",
            )
            (d / "ground.csv").write_text(
                "sensor_id,latitude_deg,longitude_deg,altitude_m,max_solar_elevation_deg\n"
                "rennes,48.1173,-1.6778,45,-6\n",
                encoding="utf-8",
            )
            out = d / "vis.csv"
            try:
                compute_visibility_metrics(
                    target_3le_path=d / "targets.3le",
                    instants_utc=["2025-03-08T00:00:00.000Z"],
                    propagator_type="sgp4",
                    output_csv_path=out,
                    ground_sensors_csv_path=d / "ground.csv",
                )
            except RuntimeError as exc:
                self.assertIn("Orekit", str(exc))
            else:
                self.assertTrue(out.exists())


if __name__ == "__main__":
    unittest.main()
