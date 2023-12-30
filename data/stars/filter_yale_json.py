

import csv
import json





def parse_ra(ra_h, ra_m, ra_s):
    hours = float(ra_h) + float(ra_m)/60 + float(ra_s)/3600
    ra = 15 * hours
    return round(ra, 4)

def parse_dec(de_sign, de_d, de_m, de_s):
    sign = -1 if de_sign == '-' else 1
    dec = sign * (float(de_d) + float(de_m)/60 + float(de_s)/3600)
    return round(dec, 4)

res = []
with open('YaleBrightStarCatalog/bsc5-all.json') as f:

    yale = json.load(f)

    for star in yale:
        ra = parse_ra(star['RAh'], star['RAm'], star['RAs'])
        dec = parse_dec(star['DE-'], star['DEd'], star['DEm'], star['DEs'])
        mag = float(star['Vmag'])
        name = star['Common'] if 'Common' in star else ''

        if mag < 6:
            #res.append({'ra':ra, 'dec':dec, 'mag':mag, 'name':name})

            out_star = [ra, dec, mag]
            if name:
                out_star.append(name)

            res.append(out_star)

print(json.dumps(res, separators=(',', ':')))

