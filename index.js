
// https://en.wikipedia.org/wiki/Julian_day
// https://github.com/mourner/suncalc/blob/master/suncalc.js 
const dayMs = 1000 * 60 * 60 * 24,
    J1970 = 2440588.0,
    J2000 = 2451545.0,
    DAYS_IN_CENTURY = 36525.0;

function julianCenturies(jd) { return (jd - J2000) / DAYS_IN_CENTURY; }

function toJulian(date) { return date.valueOf() / dayMs - 0.5 + J1970; }


function mod(m, n) {
    return ((m%n)+n)%n;
};


const rad = (deg) => deg * Math.PI / 180
const degree = (radian) => radian * (180 / Math.PI)
const sin = (deg) => Math.sin(rad(deg)), 
      cos = (deg) => Math.cos(rad(deg)), 
      tan = (deg) => Math.tan(rad(deg)), 
      acos = (x) => degree(Math.acos(x)), 
      asin = (x) => degree(Math.asin(x)), 
      atan = (x) => degree(Math.atan(x)), 
      atan2 = (x, y) => degree(Math.atan2(x, y)), 
      PI = Math.PI; 


// https://www.movable-type.co.uk/scripts/gis-faq-5.1.html
// returns angle of arc subtended by earth
// returns non-negatives value
function haversineDist(p1, p2) {
    const dlon = p2.long - p1.long
    const dlat = p2.lat - p1.lat
    const a = sin(dlat/2)**2 + cos(p1.lat) * cos(p2.lat) * sin(dlon/2)**2;
    return 2 * asin(Math.min(1, Math.sqrt(a)));
}

// since we are looking for the place at solar noon,  
// and hour angle H = 0 = side real time - right ascension, side real time == ra  
// theta(sidereal) = [theta0 + theta1 * (JD - J2000) - lw] mod 360 
// (A + B) mod C = (A mod C + B mod C) mod C 
// ra in degrees
function raToLong(jd, ra) { 
    return (280.1470 + 360.9856235 * (jd - J2000) - ra) % 360 
}

function eqCelestialToLatLong(jd, ra, dec) {
    return { lat: dec, long: raToLong(jd, ra) } 
}


// https://www.movable-type.co.uk/scripts/latlong.html - Bearing
// http://mathforum.org/library/drmath/view/55417.html
function bearing(p1, p2) {
    const y = sin(p2.long - p1.long) * cos(p2.lat);
    const x = cos(p1.lat) * sin(p2.lat) - sin(p1.lat) * cos(p2.lat) * cos(p2.long - p1.long);
    const theta = atan2(y, x)

    // since using long west, x is negative to normally long
    // this means 0->180 will be from -y axis clockwise to +y axis, and 0 -> -180 will be mapped to -y axis to +y axis ccw
    // hence to get to standard compass degrees, need to subtract 90, then do appropiate mod
    return mod(-theta, 360);
}

// computer altitude/azimuth from ra/dec of a celestial object and its parallax angle and current time, and user location in lat long
// celestial object is in equitorial celestial frame (not ecliptic)
// assume parallax is 0
function getAltAz(jd, userLoc, celestial) {
    userLoc = toLatLongWest(userLoc) 
    // celestial ra from hvg dataset is in hours, so times 15 to get degrees
    const latLongUnder = eqCelestialToLatLong(jd, celestial.ra * 15, celestial.dec)
    console.log('jd', jd)
    console.log('userLoc', userLoc)
    console.log('celestial', celestial)
    console.log('latLongUnder', latLongUnder)

    const angleDist = haversineDist(latLongUnder, userLoc)

    console.log('angleDist', angleDist)
    const altitude = 90 - angleDist
    const azimuth = bearing(userLoc, latLongUnder)
    return { altitude, azimuth } 
}

function toLatLongWest(latLong) {
    const long = latLong.long
    const longWest = long < 0 ? -long : 360 - long
    return {lat: latLong.lat, long: longWest } 
}

const austin = { lat: 30.2789, long: -97.7487 } 
const starSpica = { ra:13.419883, dec:-11.161322 } 
const starHamal = { ra:2.119555, dec:23.462423 } 
const starMirach= { ra:1.162194, dec:35.620558 } 

console.log(getAltAz(toJulian(Date.now()), austin, starMirach)) 
