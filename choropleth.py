#!/usr/bin/env python
#Originally written by Tyler Hartley ( http://beneathdata.com/how-to/visualizing-my-location-history/ )
#Ray Sharma 2015
'''
Creates a choropleth map from latitude/longitude data and binned shapefile. Each coordinate is assigned to a bin and each bin is assigned a color based on counts.
General usage:
    import choropleth as ch
    patch_count, patch_name, point_name = ch.choropleth(shapefilename, title, suptitle, savename, figwidth)

This code breaks down for small samples (<100 counts per bin) due to hardcoded rounding in create_patches (breaks = list(np.linspace(0., maximum-maximum%100, 9)) + [1e20]).
This can be fixed by setting a new modulus to round down by.
'''

import numpy as np
from shapely.geometry import Point, Polygon, MultiPoint, MultiPolygon
from shapely.prepared import prep
from matplotlib.collections import PatchCollection
from descartes import PolygonPatch
import matplotlib.pyplot as plt

def name_of_contained_points(city_point, apolygon, names):
    contained = [names[i] for i in range(len(apolygon)) if prep(apolygon[i]).contains(city_point)]
    return contained

def num_of_contained_points(apolygon, city_points):
    contained = [i for i in city_points if prep(apolygon).contains(i)]
    return int(len(contained))

def self_categorize(entry, breaks):
    for i in range(len(breaks)-1):
        if entry > breaks[i] and entry <= breaks[i+1]:
            return i
    return -1

def custom_colorbar(cmap, ncolors, breaks, **kwargs):
    from matplotlib.colors import BoundaryNorm
    from matplotlib.cm import ScalarMappable
    import matplotlib.colors as mplc

    breaklabels = ['No Counts']+["> %d counts"%(perc) for perc in breaks[:-1]]

    norm = BoundaryNorm(range(0, ncolors), cmap.N)
    mappable = ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable, **kwargs)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors+1)+0.5)
    colorbar.set_ticklabels(range(0, ncolors))
    colorbar.set_ticklabels(breaklabels)
    return colorbar

def open_shape(shapefilename): #Read in shapefile with fiona, find bounds of GIS data
    import fiona
    shp = fiona.open(shapefilename+'.shp')
    coords = shp.bounds
    return coords

def create_patches(shapefilename, latitude, longitude, coords): #Create mpl_toolkits basemap for projecting shapefile
    from mpl_toolkits.basemap import Basemap
    w, h = coords[2] - coords[0], coords[3] - coords[1]
    extra = 0.01

    m = Basemap(
        projection='tmerc', ellps='WGS84',
        lon_0=np.mean([coords[0], coords[2]]),
        lat_0=np.mean([coords[1], coords[3]]),
        llcrnrlon=coords[0] - extra * w,
        llcrnrlat=coords[1] - (extra * h),
        urcrnrlon=coords[2] + extra * w,
        urcrnrlat=coords[3] + (extra * h),
        resolution='i',  suppress_ticks=True)
    _out = m.readshapefile(shapefilename, name='map', drawbounds=False, color='none', zorder=2)

    #Construct polygons from shapefile, pull names
    poly  = [Polygon(hood_points) for hood_points in m.map]
    name = [hood['S_HOOD'] for hood in m.map_info]

    #Get location information from incoming data and match with polygon patches
    all_points = MultiPoint([Point(m(mapped_x, mapped_y)) for mapped_x, mapped_y in zip(longitude, latitude)])
    hood_polygons = prep(MultiPolygon(poly))
    city_points = [i for i in all_points if hood_polygons.contains(i)]

    #Create patches for each region or "hood"
    patches = np.array([PolygonPatch(i,  ec='#555555', lw=.8, alpha=1., zorder=4) for i in poly])
    pc = PatchCollection(patches, match_original=True)

    #Create patches with counts for each region
    hood_count = np.array([num_of_contained_points(i, city_points) for i in poly])
    point_names = np.array([name_of_contained_points(i, poly, name) for i in all_points])

    maximum = hood_count.max()
    breaks = list(np.linspace(0., maximum-maximum%100, 9)) + [1e20]

    #Create bins for labeling colorbar
    jenk_bins = np.array([self_categorize(i, breaks) for i in hood_count])

    return pc, jenk_bins,  m, breaks, hood_count, name, point_names

def set_patchcolor(pc, cmap, jenk_bins): #Create colormap & set patch color through normalized patch counts
    jenk_norm = (jenk_bins - jenk_bins.min()) / (jenk_bins.max() - jenk_bins.min())
    cmap_list = [cmap(val) for val in jenk_norm]
    pc.set_facecolor(cmap_list)
    return pc

def plotter(pc, m, coords, cmap, breaks, figwidth, title, suptitle, savename): #Create figure, plot patch collection
    w, h = coords[2] - coords[0], coords[3] - coords[1]
    extra = 0.01
    fig = plt.figure(1, figsize=(figwidth, figwidth*h/w))
    ax = fig.add_subplot(111, axisbg='#555555', frame_on=False)
    ax.add_collection(pc)

    #Draw a map scale
    m.drawmapscale(coords[0]+0.12, coords[1]+0.008,
        coords[0], coords[1], 5.,
        units='km', fontsize=16, barstyle='fancy', labelstyle='simple',
        fillcolor1='w', fillcolor2='#555555', fontcolor='#555555',
        zorder=5, ax=ax, )

    cbar  = custom_colorbar(cmap, len(breaks)+1, breaks, shrink=0.5)
    cbar.ax.tick_params(labelsize=10)

    ax.set_title(title, fontsize=10, y=0.995)
    fig.suptitle(suptitle, fontsize=15, fontweight='bold', y=0.92)
    plt.savefig(savename, dpi=100, frameon=False, bbox_inches='tight', pad_inches=0.5, facecolor='#F2F2F2')

def choropleth(shapefilename, cmap, latitude, longitude, title, suptitle, savename, figwidth): #Outputs choropleth map
    boundaries = open_shape(shapefilename)
    pc, color_bins, m, breaks, hood_count, name, point_names  = create_patches(shapefilename, latitude,
                                                    longitude, boundaries)
    pc = set_patchcolor(pc, cmap, color_bins)
    plotter(pc, m, boundaries, cmap, breaks, figwidth, title, suptitle, savename)
    return hood_count, name, point_names
