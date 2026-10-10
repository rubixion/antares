"""Builds ../data.js from three public catalogs.

  HYG v4.1 star catalog         https://github.com/astronexus/HYG-Database (CC BY-SA 4.0)
  Stellarium "modern" sky culture https://github.com/Stellarium/stellarium (figures, names, IAU boundaries)
  NASA Exoplanet Archive        https://exoplanetarchive.ipac.caltech.edu (pscomppars, every confirmed planet)

Usage:  python tools/build_data.py [cache_dir]
Downloads are cached in cache_dir (default tools/.cache). Delete it to refresh the snapshot.
"""
import csv, json, math, os, sys, urllib.parse, urllib.request
from collections import defaultdict
from datetime import date

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '.cache')
OUT = os.path.join(HERE, '..', 'data.js')
MAG_LIMIT = 6.5  # naked-eye limit under a dark sky

PLANET_COLS = ('pl_name,hostname,hip_name,sy_snum,pl_orbsmax,pl_orbper,pl_orbeccen,pl_rade,pl_bmasse,'
               'pl_bmassprov,pl_eqt,pl_insol,disc_year,discoverymethod,disc_facility,st_spectype,st_teff,'
               'st_rad,st_mass,st_lum,sy_dist,sy_vmag,ra,dec')
SOURCES = {
    'hyg.csv': 'https://raw.githubusercontent.com/astronexus/HYG-Database/main/hyg/CURRENT/hygdata_v41.csv',
    'modern.json': 'https://raw.githubusercontent.com/Stellarium/stellarium/master/skycultures/modern/index.json',
    'planets.json': 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=json&query='
                    + urllib.parse.quote(f'select {PLANET_COLS} from pscomppars'),
}

GREEK = dict(Alp='α', Bet='β', Gam='γ', Del='δ', Eps='ε', Zet='ζ', Eta='η', The='θ', Iot='ι', Kap='κ', Lam='λ',
             Mu='μ', Nu='ν', Xi='ξ', Omi='ο', Pi='π', Rho='ρ', Sig='σ', Tau='τ', Ups='υ', Phi='φ', Chi='χ',
             Psi='ψ', Ome='ω')
SUP = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')

GENITIVE = dict(
    And='Andromedae', Ant='Antliae', Aps='Apodis', Aqr='Aquarii', Aql='Aquilae', Ara='Arae', Ari='Arietis',
    Aur='Aurigae', Boo='Boötis', Cae='Caeli', Cam='Camelopardalis', Cnc='Cancri', CVn='Canum Venaticorum',
    CMa='Canis Majoris', CMi='Canis Minoris', Cap='Capricorni', Car='Carinae', Cas='Cassiopeiae', Cen='Centauri',
    Cep='Cephei', Cet='Ceti', Cha='Chamaeleontis', Cir='Circini', Col='Columbae', Com='Comae Berenices',
    CrA='Coronae Australis', CrB='Coronae Borealis', Crv='Corvi', Crt='Crateris', Cru='Crucis', Cyg='Cygni',
    Del='Delphini', Dor='Doradus', Dra='Draconis', Equ='Equulei', Eri='Eridani', For='Fornacis', Gem='Geminorum',
    Gru='Gruis', Her='Herculis', Hor='Horologii', Hya='Hydrae', Hyi='Hydri', Ind='Indi', Lac='Lacertae',
    Leo='Leonis', LMi='Leonis Minoris', Lep='Leporis', Lib='Librae', Lup='Lupi', Lyn='Lyncis', Lyr='Lyrae',
    Men='Mensae', Mic='Microscopii', Mon='Monocerotis', Mus='Muscae', Nor='Normae', Oct='Octantis',
    Oph='Ophiuchi', Ori='Orionis', Pav='Pavonis', Peg='Pegasi', Per='Persei', Phe='Phoenicis', Pic='Pictoris',
    Psc='Piscium', PsA='Piscis Austrini', Pup='Puppis', Pyx='Pyxidis', Ret='Reticuli', Sge='Sagittae',
    Sgr='Sagittarii', Sco='Scorpii', Scl='Sculptoris', Sct='Scuti', Ser='Serpentis', Sex='Sextantis',
    Tau='Tauri', Tel='Telescopii', Tri='Trianguli', TrA='Trianguli Australis', Tuc='Tucanae',
    UMa='Ursae Majoris', UMi='Ursae Minoris', Vel='Velorum', Vir='Virginis', Vol='Volantis', Vul='Vulpeculae')
ZODIAC = {'Ari', 'Tau', 'Gem', 'Cnc', 'Leo', 'Vir', 'Lib', 'Sco', 'Sgr', 'Cap', 'Aqr', 'Psc'}


def fetch(name):
    os.makedirs(CACHE, exist_ok=True)
    path = os.path.join(CACHE, name)
    if not os.path.exists(path):
        print('downloading', name)
        urllib.request.urlretrieve(SOURCES[name], path)
    return path


def num(s):
    try:
        v = float(s)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def rnd(v, d):
    return None if v is None else round(v, d)


def sig(v, n=3):
    """Round to n significant figures (luminosities span 1e-4..1e6)."""
    if not v:
        return v
    return float(f'{v:.{n}g}')


# Astrophysics

# Spectral-class anchors, Teff in K (Pecaut & Mamajek 2013, dwarfs). Used for O/B stars, whose B-V is
# saturated and reddened, and wherever B-V is missing.
SPEC_ANCHORS = [(3, 44900), (5, 41400), (7, 36500), (9, 31500), (10, 31400), (11, 26000), (12, 20600),
                (13, 17000), (15, 15700), (18, 12500), (19, 10700), (20, 9700), (22, 8840), (25, 8080),
                (30, 7220), (35, 6510), (40, 5920), (42, 5770), (45, 5660), (50, 5280), (52, 4990),
                (55, 4450), (60, 3850), (62, 3560), (65, 3060), (68, 2650), (69, 2400)]


def spec_code(spect):
    if not spect:
        return None
    s = spect.strip()
    while s and s[0] in 'dDsgc':  # luminosity prefixes such as "dM2", "gK0"
        s = s[1:]
    if not s or s[0] not in 'OBAFGKM':
        return None
    code = 'OBAFGKM'.index(s[0]) * 10
    digits = ''
    for ch in s[1:]:
        if ch.isdigit() or (ch == '.' and digits and '.' not in digits):
            digits += ch
        else:
            break
    return code + (float(digits.rstrip('.')) if digits else 5)


def teff_from_spec(code):
    xs, ys = zip(*SPEC_ANCHORS)
    return float(np.interp(code, xs, ys))


def teff_from_bv(bv):
    """Ballesteros (2012)."""
    bv = min(max(bv, -0.3), 2.0)
    return 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62))


def teff(spect, bv):
    code = spec_code(spect)
    if code is not None and (code < 20 or bv is None):
        return teff_from_spec(code)
    if bv is not None:
        return teff_from_bv(bv)
    return None


# Torres (2010) correction of Flower (1996): BC_V as a polynomial in log10 Teff.
_BC = [
    (3.70, [-0.190537291496456e+05, 0.155144866764412e+05, -0.421278819301717e+04, 0.381476328422343e+03]),
    (3.90, [-0.370510203809015e+05, 0.385672629965804e+05, -0.150651486316025e+05, 0.261724637119416e+04,
            -0.170623810323864e+03]),
    (9.99, [-0.118115450538963e+06, 0.137145973583929e+06, -0.636233812100225e+05, 0.147412923562646e+05,
            -0.170587278406872e+04, 0.788731721804990e+02]),
]


def bolometric_correction(t):
    x = math.log10(t)
    coeffs = next(c for limit, c in _BC if x < limit)
    return sum(c * x ** i for i, c in enumerate(coeffs))


def lum_bol(absmag, t):
    return 10 ** (-0.4 * (absmag + bolometric_correction(t) - 4.74))


def radius(l, t):
    return math.sqrt(l) * (5772 / t) ** 2


def radius_from_mass(m):
    """Chen & Kipping (2017) forecaster, mean relation, Earth units."""
    if m < 2.04:
        return m ** 0.279
    if m < 132:
        return 0.808 * m ** 0.589
    return 17.7 * m ** -0.044


# Coordinates

def vec(ra_deg, dec_deg):
    a, d = math.radians(ra_deg), math.radians(dec_deg)
    return np.array([math.cos(d) * math.cos(a), math.cos(d) * math.sin(a), math.sin(d)])


def radec(v):
    v = v / np.linalg.norm(v)
    return math.degrees(math.atan2(v[1], v[0])) % 360, math.degrees(math.asin(max(-1, min(1, v[2]))))


def precession_matrix_j2000_to(jd):
    """IAU 1976 (Lieske) precession, J2000.0 mean equator to the mean equator of `jd`."""
    t = (jd - 2451545.0) / 36525
    sec = math.radians(1 / 3600)
    zeta = (2306.2181 * t + 0.30188 * t ** 2 + 0.017998 * t ** 3) * sec
    z = (2306.2181 * t + 1.09468 * t ** 2 + 0.018203 * t ** 3) * sec
    theta = (2004.3109 * t - 0.42665 * t ** 2 - 0.041833 * t ** 3) * sec

    def r3(p):
        c, s = math.cos(p), math.sin(p)
        return np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    def r2(p):
        c, s = math.cos(p), math.sin(p)
        return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])

    return r3(-z) @ r2(theta) @ r3(-zeta)


P_TO_1875 = precession_matrix_j2000_to(2405889.258550475)  # B1875.0
P_FROM_1875 = P_TO_1875.T


def hms(s):
    h, m, sec = (float(x) for x in s.split(':'))
    return (h + m / 60 + sec / 3600) * 15


def dms(s):
    sign = -1 if s.startswith('-') else 1
    d, m, sec = (float(x) for x in s.lstrip('+-').split(':'))
    return sign * (d + m / 60 + sec / 3600)


def parse_edges(edges):
    """IAU boundary segments in the B1875 frame: (type, ra1, dec1, ra2, dec2, conA, conB)."""
    out = []
    for e in edges:
        p = e.split()
        out.append((p[1][0], hms(p[2]), dms(p[3]), hms(p[4]), dms(p[5]),
                    'Ser' if p[6].startswith('SER') else p[6], 'Ser' if p[7].startswith('SER') else p[7]))
    return out


class Boundaries:
    """Point-in-constellation for J2000 positions, using the IAU B1875 boundaries.

    Boundaries run along B1875 meridians and parallels, so a ray cast south along a meridian only
    crosses parallel ('P') segments. The ray ends at the south celestial pole, inside Octans; each
    crossing flips membership of the two constellations sharing that edge, and the one left
    flipped on is the constellation that holds the point.
    """

    def __init__(self, edges, abbrs):
        self.idx = {a.upper(): i for i, a in enumerate(abbrs)}
        self.n = len(abbrs)
        self.oct = self.idx['OCT']
        pe = [e for e in edges if e[0] == 'P']
        self.dec = np.array([e[2] for e in pe])
        lo, span = [], []
        for _, ra1, _, ra2, _, _, _ in pe:
            d = (ra2 - ra1) % 360
            if d > 180:  # take the short way round
                ra1, d = ra2, 360 - d
            lo.append(ra1)
            span.append(d)
        self.lo, self.span = np.array(lo), np.array(span)
        self.a = np.array([self.idx[e[5].upper()] for e in pe])
        self.b = np.array([self.idx[e[6].upper()] for e in pe])

    def lookup(self, ra_j2000, dec_j2000):
        """Vectorised: arrays of J2000 degrees in, constellation indices out."""
        v = np.stack([np.cos(np.radians(dec_j2000)) * np.cos(np.radians(ra_j2000)),
                      np.cos(np.radians(dec_j2000)) * np.sin(np.radians(ra_j2000)),
                      np.sin(np.radians(dec_j2000))])
        w = P_TO_1875 @ v
        ra = np.degrees(np.arctan2(w[1], w[0])) % 360
        dec = np.degrees(np.arcsin(np.clip(w[2], -1, 1)))
        parity = np.zeros((self.n, ra.size), dtype=bool)
        parity[self.oct] = True
        for k in range(self.dec.size):
            hit = (self.dec[k] < dec) & (((ra - self.lo[k]) % 360) < self.span[k])
            parity[self.a[k]] ^= hit
            parity[self.b[k]] ^= hit
        assert (parity.sum(axis=0) == 1).all(), 'boundary parity broken'
        return parity.argmax(axis=0)


def boundary_polylines(edges):
    """Precessed J2000 polylines. Parallels are small circles, so interpolate them; meridians stay great
    circles under precession, and a straight chord seen from the sphere's centre already draws one."""
    out = []
    for kind, ra1, dec1, ra2, dec2, a, b in edges:
        if kind == 'P':
            d = (ra2 - ra1) % 360
            if d > 180:
                d -= 360
            n = max(2, int(abs(d) / 2) + 2)
            pts = [(ra1 + d * i / (n - 1), dec1) for i in range(n)]
        else:
            pts = [(ra1, dec1), (ra2, dec2)]
        flat = []
        for ra, dec in pts:
            r, d_ = radec(P_FROM_1875 @ vec(ra, dec))
            flat += [round(r, 2), round(d_, 2)]
        out.append([a, b, flat])
    return out


# Build

def main():
    hyg = list(csv.DictReader(open(fetch('hyg.csv'), encoding='utf8')))
    sky = json.load(open(fetch('modern.json'), encoding='utf8'))
    planets = json.load(open(fetch('planets.json'), encoding='utf8'))

    cons_src = sky['constellations']
    abbrs = [c['id'].split()[-1] for c in cons_src]
    assert len(abbrs) == 88 and set(abbrs) == set(GENITIVE), 'expected the 88 IAU constellations'
    con_index = {a: i for i, a in enumerate(abbrs)}
    line_hips = {h for c in cons_src for line in c['lines'] for h in line}
    names = {int(k.split()[1]): v[0]['english'] for k, v in sky['common_names'].items() if k.startswith('HIP ')}

    edges = parse_edges(sky['edges'])
    bounds = Boundaries(edges, abbrs)

    # Exoplanet systems, grouped by host
    methods, facilities = [], []

    def intern(table, s):
        if s not in table:
            table.append(s)
        return table.index(s)

    by_host = defaultdict(list)
    for p in planets:
        by_host[p['hostname']].append(p)

    systems = []
    host_by_hip = {}
    for host, ps in sorted(by_host.items()):
        s0 = ps[0]
        st_mass = num(s0['st_mass'])
        st_teff = num(s0['st_teff'])
        st_rad = num(s0['st_rad'])
        log_l = num(s0['st_lum'])
        lum = 10 ** log_l if log_l is not None else (st_rad ** 2 * (st_teff / 5772) ** 4 if st_rad and st_teff else None)
        rows = []
        for p in sorted(ps, key=lambda p: (num(p['pl_orbsmax']) or num(p['pl_orbper']) or 1e9)):
            a, per = num(p['pl_orbsmax']), num(p['pl_orbper'])
            r, m, teq = num(p['pl_rade']), num(p['pl_bmasse']), num(p['pl_eqt'])
            flags, mstar = 0, st_mass or 1.0
            if a is None and per is not None:
                a, flags = (mstar * (per / 365.25) ** 2) ** (1 / 3), flags | 1
            if per is None and a is not None:
                per, flags = 365.25 * math.sqrt(a ** 3 / mstar), flags | 2
            if a is None:
                continue  # no orbit at all (a handful of microlensing/timing entries)
            if r is None and m is not None:
                r, flags = radius_from_mass(m), flags | 4
            if teq is None and lum:
                teq, flags = 278.6 * lum ** 0.25 / math.sqrt(a), flags | 8
            rows.append([p['pl_name'], sig(a, 4), sig(per, 5), rnd(num(p['pl_orbeccen']), 3), rnd(r, 2),
                         sig(m, 4), 1 if p['pl_bmassprov'] == 'Msini' else 0, rnd(teq, 0),
                         sig(num(p['pl_insol']), 3), p['disc_year'],
                         intern(methods, p['discoverymethod']), intern(facilities, p['disc_facility']), flags])
        if not rows:
            continue
        hip = int(s0['hip_name'].split()[1]) if s0['hip_name'] else 0
        systems.append(dict(host=host, hip=hip, ra=num(s0['ra']), dec=num(s0['dec']), dist=num(s0['sy_dist']),
                            vmag=num(s0['sy_vmag']), spect=s0['st_spectype'] or '', teff=st_teff, rad=st_rad,
                            mass=st_mass, lum=lum, snum=s0['sy_snum'], planets=rows))
        if hip and hip not in host_by_hip:
            host_by_hip[hip] = len(systems) - 1

    sys_con = bounds.lookup(np.array([s['ra'] for s in systems]), np.array([s['dec'] for s in systems]))

    # Stars
    stars, star_by_hip = [], {}
    for r in hyg:
        mag, dist = num(r['mag']), num(r['dist'])
        hip = int(r['hip']) if r['hip'] else 0
        if r['proper'] == 'Sol' or mag is None or not dist:
            continue
        if mag > MAG_LIMIT and hip not in line_hips:
            continue
        ci, absmag = num(r['ci']), num(r['absmag'])
        t = teff(r['spect'], ci) or 5800
        known = dist < 100000  # HYG uses 100000 pc for "no parallax"
        l = lum_bol(absmag, t) if known and absmag is not None else None
        bayer = r['bayer']
        greek = ''
        if bayer:
            letter, _, n = bayer.partition('-')
            greek = GREEK.get(letter, letter) + n.translate(SUP)
        con = r['con']
        if con not in con_index:
            continue
        sys_i = host_by_hip.get(hip, -1)
        if sys_i >= 0:  # prefer the archive's measured stellar parameters for planet hosts
            s = systems[sys_i]
            if s['teff']:
                t = s['teff']
            if s['lum']:
                l = s['lum']
        rad = (systems[sys_i]['rad'] if sys_i >= 0 and systems[sys_i]['rad'] else (radius(l, t) if l else None))
        star_by_hip[hip] = len(stars)
        stars.append([hip, round(float(r['ra']) * 15, 4), round(float(r['dec']), 4), mag, rnd(ci, 2),
                      rnd(dist, 2) if known else None, r['spect'], names.get(hip) or r['proper'] or '',
                      greek, r['flam'] or '', con_index[con], int(r['hd']) if r['hd'] else 0,
                      int(r['hr']) if r['hr'] else 0, round(t), sig(l, 3), sig(rad, 3), sys_i])
    for s in systems:
        s['star'] = star_by_hip.get(s['hip'], -1) if s['hip'] and host_by_hip.get(s['hip']) is not None else -1
    for i, s in enumerate(systems):
        if s['star'] >= 0 and stars[s['star']][16] != i:
            s['star'] = -1  # second host sharing a HIP number (binary components) gets its own marker

    # Constellations
    cons = []
    star_count = defaultdict(int)
    for st in stars:
        if st[3] <= MAG_LIMIT:
            star_count[st[10]] += 1
    for ci_, c in enumerate(cons_src):
        a = abbrs[ci_]
        lines = [[star_by_hip[h] for h in line] for line in c['lines']]
        figure = {i for line in lines for i in line}
        anchor = radec(sum(vec(stars[i][1], stars[i][2]) for i in figure))
        members = [i for i, st in enumerate(stars) if st[10] == ci_]
        brightest = min(members, key=lambda i: stars[i][3])
        cons.append(dict(abbr=a, name=c['common_name']['native'], meaning=c['common_name']['english'],
                         gen=GENITIVE[a], zodiac=a in ZODIAC, lines=lines, ra=round(anchor[0], 2),
                         dec=round(anchor[1], 2), bright=brightest, count=star_count[ci_],
                         hosts=int((sys_con == ci_).sum())))

    # 1 degree lookup grid: which constellation each direction belongs to (for clicks and the HUD)
    gra, gdec = np.meshgrid(np.arange(360) + 0.5, np.arange(180) - 89.5)
    grid = bounds.lookup(gra.ravel(), gdec.ravel())
    alphabet = ''.join(chr(c) for c in range(35, 127) if c != 92)[:88]
    grid_str = ''.join(alphabet[i] for i in grid)

    out_systems = [[s['host'], s['hip'], round(s['ra'], 4), round(s['dec'], 4), rnd(s['dist'], 2),
                    rnd(s['vmag'], 2), s['spect'], s['teff'], s['rad'], s['mass'], sig(s['lum'], 3), s['snum'],
                    int(sys_con[i]), s['star'], s['planets']] for i, s in enumerate(systems)]

    data = dict(
        built=date.today().isoformat(),
        stars=stars, cons=cons, bounds=boundary_polylines(edges), grid=grid_str, alphabet=alphabet,
        systems=out_systems, methods=methods, facilities=facilities,
    )
    js = ('// Generated by tools/build_data.py. Do not edit by hand.\n'
          '// Stars: HYG v4.1 (CC BY-SA 4.0). Figures, names and IAU boundaries: Stellarium (GPL-2.0).\n'
          '// Planets: NASA Exoplanet Archive, Planetary Systems Composite Parameters.\n'
          'window.SKY=' + json.dumps(data, ensure_ascii=False, separators=(',', ':')) + ';\n')
    with open(OUT, 'w', encoding='utf8') as f:
        f.write(js)

    # Self-checks
    hyg_con = np.array([st[10] for st in stars])
    got = bounds.lookup(np.array([st[1] for st in stars]), np.array([st[2] for st in stars]))
    agree = (got == hyg_con).mean()
    print(f'stars {len(stars)}  systems {len(systems)}  planets {sum(len(s["planets"]) for s in systems)}'
          f'  boundary agreement with HYG {agree:.4%}  size {len(js) / 1e6:.2f} MB')
    assert agree > 0.995, 'boundary lookup disagrees with HYG constellation labels'


def self_test():
    # 1875 (0h, 0 deg) precesses to roughly (0h06m24s, +0d41.7m) at J2000
    ra, dec = radec(P_FROM_1875 @ vec(0, 0))
    assert abs(ra - 1.60) < 0.02 and abs(dec - 0.695) < 0.01, (ra, dec)
    assert abs(bolometric_correction(5772) + 0.07) < 0.06, bolometric_correction(5772)
    assert abs(lum_bol(4.83, 5772) - 1) < 0.08
    assert abs(teff_from_bv(0.65) - 5778) < 30
    assert abs(teff('B1.5IV+...', -0.23) - 23300) < 1
    assert abs(radius_from_mass(1) - 1) < 1e-9 and 10 < radius_from_mass(318) < 15


if __name__ == '__main__':
    self_test()
    main()
