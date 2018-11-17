# /usr/bin/env python
import urllib
import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
from matplotlib.colors import LinearSegmentedColormap
import time
import datetime
from math import sin, cos, atan2, radians
import os


def f_prev(url):
    all = urllib.request.urlopen(url).read().decode()
    coor = all.split('"coords":{')[1].split('}}]\'>')[0].replace(',"lon"','').split(':')[1:]
    coor = [float(c) for c in coor]
    all = all.split('forecast: ')[1].split(',\n\t\t\t')[0]
    all = eval(all.replace('false','\'false\'').replace('true','\'true\''))
    prev_vague = pd.concat([pd.DataFrame({all[i]['timestamp'] : all[i]['swell']}).T for i in range(len(all))],sort=True); prev_vague.columns = ['S_'+c for c in prev_vague.columns]; prev_vague.drop(['S_components'],axis=1)
    prev_vent = pd.concat([pd.DataFrame({all[i]['timestamp'] : all[i]['wind']}).T for i in range(len(all))],sort=True); prev_vent.columns = ['W_'+c for c in prev_vent.columns]
    prev_rate = pd.concat([pd.DataFrame({all[i]['timestamp'] : all[i]}).T for i in range(len(all))],sort=True); prev_rate.columns = ['R_'+c for c in prev_rate.columns]
    prev = pd.concat([prev_vague,prev_vent,prev_rate],axis=1)
    return prev,coor

def f_loc(spot):
    api_key = 'AIzaSyAjRoFrgc0_qLLi8XEy8QVHIkq8otnaEpk'
    r = requests.get('https://maps.googleapis.com/maps/api/geocode/json?address='+spot.replace('-','+')+'+surf+France&key='+api_key); r = r.json()
    loc = list(r['results'][0]['geometry']['location'].values())
    return loc

def f_add():
	url = 'https://magicseaweed.com/site-map.php'
	all = str(urllib.request.urlopen(url).read()).split('optgroup label="France">')[1]
	add = ['https://magicseaweed.com'+i[0:i.find('"')] for i in all[0:all.find('</optgroup')].split('value="')[1:]]
	spots = [addi[25:addi.find('-Surf-Report')]  for addi in add]
	return dict(zip(spots,add))

def f_get_all(spots=[]):
    add = f_add()
    if spots == [] : spots = list(add.keys())
    loc = dict(zip(spots,[[]]*len(spots)))
    prev = dict(zip(spots,[[]]*len(spots)))
    erreurs = []; succes =[] ; t = 0
    force = ['Siouville','Les-Sables-dOlonne','Les-Conches-Bud-Bud','Biarritz-Grande-Plage','Etretat','La-Torche']
    while (len([e for e in force if e in erreurs]) or succes == []) and t <= 2:
        #if t > 0 : print('Nouvelle passe...')
        for spot in spots:
            try:
                #loc[spot] = f_loc(spot)
                prev[spot],loc[spot] = f_prev(add[spot])
                succes += [spot]
                if t > 0 : erreurs.remove(spot)
            except:
                if t == 0:  erreurs += [spot]
            print('\r',end='\r')
            print(len(succes),' succes et ',len(erreurs),' echec(s)',end='\r')
        for s in succes :
            if s in spots : spots.remove(s)
        t += 1
    for e in erreurs :
        prev.pop(e)
        loc.pop(e)
    loc = pd.DataFrame.from_dict(loc,orient='index')
    prev = pd.concat(list(prev.values()),keys=list(prev.keys()),names=['Spot'],axis=1).swaplevel(0,1,1)
    return add,loc,prev,erreurs

def f_lit(loc,path):
    coor = loc.copy()
    file = os.path.join(path,'litoral.csv')
    lit = pd.read_csv(file,sep=';'); lit['dist'] = np.cumsum([0]+[np.sqrt((lit['lat'][i-1]-lit['lat'][i])**2 + (lit['long'][i-1]-lit['long'][i])**2) for i in list(lit.index)[1:]]); lit['dist'] = lit['dist']/max(lit['dist'])       
    lit_all = lit.set_index(['dist']); lit_all = lit_all.reindex(sorted(list(set(list(lit_all.index)+list(np.arange(0,1,1e-4)))))).interpolate('index')
    spot_s = loc.copy()
    spot_s['x'] = [(np.sqrt((lit_all['lat']-lat)**2+(lit_all['long']-long)**2)).idxmin() for lat,long in np.array(coor)]
    spot_s['x1'] = spot_s['x']
    spot_s['d'] = np.array([min(np.sqrt((lit['lat']-lat)**2+(lit['long']-long)**2)) for lat,long in np.array(coor)])*60/1.85
    spot_f = spot_s[spot_s['d']<5].reset_index().set_index('x').sort_index()
    spot_f['d_cumul'] = spot_f['x1'] #list(np.cumsum(spot_f['d']))-spot_f['d'].iloc[0]
    return spot_f

def f_post(spot_f,prev,path,Max=1e6):
    plt.close('all')
    class MidpointNormalize(colors.Normalize):
        def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
            self.midpoint = midpoint
            super().__init__(vmin, vmax, clip)
        def __call__(self, value, clip=None):
            x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
            return np.ma.masked_array(np.interp(value, x, y))
    cmap=LinearSegmentedColormap.from_list('rg',list(zip(list(np.linspace(0,1,6)),["darkred","red","yellow","yellow","green","darkgreen"])), N=256)
    cmap_r=LinearSegmentedColormap.from_list('gr',list(zip(list(np.linspace(0,1,6)),["darkred","red","yellow","yellow","green","darkgreen"][::-1])), N=256)
    data = [[]]*5;
    spot_f = spot_f[spot_f['d_cumul']<Max]
    top2bot = [s for s in list(spot_f['index']) if s in prev.columns.levels[1]]
    titres = ['Note', 'Vagues (m)','Periode (s)','Vent (km/h)','Vent offshore (%)']
    prev = prev.replace('false',0)
    data[1] = np.array(prev['S_absMinBreakingHeight'][top2bot]); #cmaps[0]='jet_r';vmaxs[0]
    data[2] = np.array(prev['S_period'][top2bot])
    data[3] = np.array(prev['W_speed'][top2bot])
    data[4] = abs(np.array(prev['W_direction'][top2bot]-prev['S_direction'][top2bot]))/3.6
    data[0] = np.array(prev['R_solidRating'][top2bot])
    dist = list(spot_f['d_cumul'])
    X_ticks = [spot_f.set_index('index').loc[s,'d_cumul'] for s in [spot_f['index'].iloc[0],'Etretat','Siouville','La-Torche','Trestraou','Les-Sables-dOlonne',spot_f['index'].iloc[-1]] if s in list(spot_f['index'])]
    T_ticks = list(prev.index)[::8]+[list(prev.index)[-1]]
    t,x = np.meshgrid(np.array(prev.index),np.array(dist))
    fig, grid = plt.subplots(nrows=5, ncols=1, sharex=False, sharey=True,figsize=(8,15))
    args = [{'vmin':0,'vmax':5,'cmap':cmap},
            {'vmin':0,'vmax':1.5,'cmap':cmap},
            {'vmin':0,'vmax':15,'cmap':cmap},
            {'vmin':0,'vmax':35,'cmap':cmap_r},
            {'vmin':0,'vmax':100,'cmap':cmap}]
    f = ['%.1f', '%.1f','%.1f','%d', '%d']
    for i,ax in enumerate(grid.flat):
        ticks = [np.linspace(round(np.array(data[i]).min(),1),round(np.array(data[i]).max(),1),7),
                 np.linspace(round(np.array(data[i]).min(),1),round(np.array(data[i]).max(),1),7),
                 np.linspace(round(np.array(data[i]).min(),0),round(np.array(data[i]).max(),0),7),
                 np.linspace(round(np.array(data[i]).min(),0),round(np.array(data[i]).max(),0),7),
                 np.linspace(0,100,7)]
        ax.set_title(titres[i])
        cb = ax.contourf(x,t,data[i].T,50,**args[i])
        ax.set_xticks(X_ticks);ax.set_yticks(T_ticks);ax.set_xticks(list(spot_f['d_cumul']), minor=True)
        ax.xaxis.grid(True, zorder=0, color ='k'); ax.yaxis.grid(True, zorder=0, color ='k')
        ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: time.strftime('%a %d, %Hh',time.localtime(x))))
        ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda x, p: top2bot[list(dist).index(x)]))
        plt.sca(ax); plt.xticks(rotation=15,ha='right'); plt.colorbar(cb,ax=ax,format=f[i],ticks=ticks[i])
    plt.suptitle('report de : '+str(datetime.datetime.now()).split(' ')[1].split('.')[0])
    fig.subplots_adjust(left=0.20, bottom=0.06, right=1.0, top=0.95, wspace=0.0, hspace=0.5)
    img_path= os.path.join(path,'report.png')
    fig.savefig(img_path, dpi=100)
    plt.close('all')
    return prev,spot_f

def f_best(prev,spots,path,ville='Paris-France'):
    #recu = recu.split(' ')
    #ville = recu[-2]; R = int(recu[-3]); nb = int(recu[-4]); j = recu[-1].split('-')
    #spots = pd.read_csv('spots.csv').set_index('index')
    #prev = pd.read_csv('notes.csv').set_index('Unnamed: 0').sort_index()
    prev = prev['R_solidRating'].sort_index()
    spots = spots.reset_index().set_index('index')
    latlong = f_loc(ville)[::-1]
    lat1 = radians(latlong[0]); lon1 = radians(latlong[1])
    spots['dlon'] = spots[1].apply(radians)-lon1
    spots['dlat'] = spots[0].apply(radians)-lat1
    spots['a'] = (((spots['dlat']/2).apply(sin))**2 +
                 (spots[0].apply(radians))*cos(lat1)*
                 ((spots['dlon']/2).apply(sin))**2)**0.5
    spots['b'] = (1-spots['a'])**0.5
    spots['dist'] = 6373.0*2*(spots['a'].apply(lambda x : atan2(x,spots[spots['a']==x]['b'])))
    prev['Spot'] = [time.strftime('%a',time.localtime(x)) for x in list(prev.index)]
    #spots_in = list(spots[spots['dist']<R].reset_index()['index'])
    spots_in = prev.set_index('Spot').columns
    table = prev.set_index('Spot')[spots_in].mean(level=0).T
    table['(km)'] = spots['dist'].apply(round)
    table = table.sort_values('(km)')
    table = table[(table.iloc[:,:-1].T >= 2).any()]
    del table.index.name
    table.columns.name = str(datetime.datetime.now()).split(' ')[1].split('.')[0]
    os.chdir(path)
    table.to_html('table.html')
    if ville == 'Paris-France' : suf = ''
    else : suf = '_'+ville.split('-')[0]
    os.system('xvfb-run wkhtmltoimage -f png --width 0 table.html table'+suf+'.png')
    os.system('rm table.html')
    return table

def push(path = os.path.join('/','home','pi','Bureau','Mirmoc_is_back')):
    os.chdir(path)
    run=['git add *.png','git commit -m "maj"','git push']
    for c in run: os.system(c)
    
def pull(path = os.path.join('/','home','pi','Bureau','Mirmoc_is_back')):
    os.chdir(path)
    run=['git fetch --all','git reset --hard origin/master','git pull origin master']
    for c in run: os.system(c)

pull()    
rep = os.path.join('/','home','pi','Bureau','Mirmoc')
Git = os.path.join('/','home','pi','Bureau','Mirmoc_is_back')
spots = ['Calais', 'Cap-Gris-Nez', 'Wimereux', 'Yport', 'Vaucottes', 'Etretat', 'Sainte-Adresse', 'Trouville', 'LAnse-du-Brick', 'Collignon', 'Siouville', 'Dielette', 'Le-Rozel', 'Surtainville', 'Hatainville', 'Plage-du-Sillon-St-Malo', 'Les-Longchamps', 'Cap-Frehel', 'Trestraou', 'Pors-Ar-Villec-Locquirec', 'Le-Dossen', 'La-Mauvaise-Greve', 'St-Pabu', 'Lampaul-Ploudalmezeau', 'Penfoul', 'Le-Gouerou', 'Blancs-Sablons', 'Porsmilin', 'Dalbosc', 'Le-Petit-Minou', 'Anse-de-Pen-hat', 'Kerloch', 'Pointe-de-Dinan', 'La-Palue', 'Cap-de-la-Chevre', 'Pors-ar-Vag', 'Les-Roches-Blanches-Pointe-Leyde', 'Porz-Theolen', 'Baie-des-Trepasses-Lescoff', 'Saint-Tugen', 'Pointe-de-Lervily', 'La-Gamelle', 'Gwendrez', 'Fouesnou', 'La-Torche', 'Porzcarn', 'Lesconil', 'Le-Kerou', 'Plage-du-Loch', 'Guidel-Les-Kaolins', 'Toulhars', 'Gavres', 'Etel', 'Penthievre', 'La-Courance', 'LErmitage', 'Gohaud', 'Noirmoutier', 'St-jean-de-monts', 'Les-Dunes', 'Sauveterre', 'LAubraie', 'Les-Sables-dOlonne', 'Saint-Nicolas', 'Les-Conches-Bud-Bud', 'Le-Lizay', 'Les-Boulassiers', 'St-Denis', 'Les-Huttes', 'Chassiron', 'Les-Allassins', 'St-Trojan', 'La-Cote-Sauvage', 'Le-Verdon', 'Pontaillac', 'Soulac', 'LAmelie', 'Le-Gurp', 'Montalivet', 'Le-Truc-Vert', 'LHorizon', 'La-Pointe', 'La-Salie', 'Messanges', 'Vieux-Boucau', 'Casernes', 'Le-Penon', 'Les-Bourdaines', 'Les-Estagnots', 'Les-Culs-Nus', 'Hossegor-La-Graviere', 'Hossegor-La-Nord', 'Hossegor-La-Sud', 'LEstacade', 'Le-Prevent', 'Le-Santocha', 'Capbreton-La-Piste-VVF', 'Labenne-Ocean', 'Ondres-Plage', 'Tarnos','La-Barre', 'Les-Cavaliers', 'La-Madrague', 'Les-Corsaires', 'Marinella', 'Sables-dOr', 'Le-Club', 'Le-VVF', 'Biarritz-Grande-Plage', 'Cote-des-Basques', 'Ilbarritz', 'Erretegia', 'Bidart', 'Parlementia', 'Guethary-Avalanche', 'Les-Alcyons', 'Guethary', 'Lafitenia', 'Erromardie', 'Sainte-Barbe', 'Ciboure-Socoa', 'Belharra-Perdun', 'Hendaye-Plage']
add,loc,prev,erreurs = f_get_all(spots)
spots = f_lit(loc,Git)
prev,spots = f_post(spots,prev,Git)
f_best(prev,spots,Git)
f_best(prev,spots,Git,'Rennes-France')
push()
