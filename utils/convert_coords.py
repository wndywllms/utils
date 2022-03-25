import astropy.units as u
import astropy.coordinates as ac

import argparse


parser = argparse.ArgumentParser(description='Covert coordinates')
parser.add_argument('lat',type=str)
parser.add_argument('lon',type=str)
parser.add_argument('-u1','--latunit',default='hourangle')
parser.add_argument('-u2','--lonunit',default='deg')

args = parser.parse_args()

lat = args.lat
lon = args.lon


c = ac.SkyCoord(lat,lon,unit=(args.latunit, args.lonunit))


print(f'{c.ra.value:.5f} {c.dec.value:.5f} degrees')
print(f'{c.ra:.5f} {c.dec:.5f}')

rahmsstr = c.ra.to_string(u.hour, precision=2)
decdmsstr = c.dec.to_string(u.degree, precision=3, alwayssign=True)
print(f'{rahmsstr} {decdmsstr}')

rahmsstr = c.ra.to_string(u.hour,sep=':', precision=2)
decdmsstr = c.dec.to_string(u.degree,sep=':', precision=3, alwayssign=True)
print(f'{rahmsstr} {decdmsstr}')

rahmsstr = c.ra.to_string(u.hour,sep=':', precision=2)
decdmsstr = c.dec.to_string(u.degree,sep='.', precision=3, alwayssign=True)
print(f'{rahmsstr} {decdmsstr}')
