

// https://en.wikipedia.org/wiki/Julian_day
// https://github.com/mourner/suncalc/blob/master/suncalc.js 
const millisPerDay = 1000 * 60 * 60 * 24
const daysPerCentury = 36525
const j2000Date = Date.UTC(2000, 0, 1, 12, 0, 0)
const j2000jd = 2451545
const toJulianCenturies = (jd) =>  (jd - j2000jd) / daysPerCentury
const toJd = (date) => (date - j2000Date) / millisPerDay + j2000jd
const mod = (m, n) => ((m%n)+n)%n 
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






const EARTH_RADIUS_KM = 6378.14

// Originally from most recent Astronomical Almanac, at least 1999-2015
// An Alternative Lunar Ephemeris Model
// https://caps.gsfc.nasa.gov/simpson/pubs/slunar.pdf
export class Moon {

    // Main public method for 
    static eclipLatLong(jd) {
        const t = toJulianCenturies(jd)
        const eclipLongBase = Moon.eclipticLongBase(t)
        const eclipLatBase = Moon.eclipticLatBase(t)
        //const [eclipLat, eclipLong] = Moon.precessionCorrection(t, eclipLatBase, eclipLongBase)
        return {
            long: eclipLongBase,
            lat: eclipLatBase
        }
    }

    static horizontalParallaxJd(jd) {
        const t = toJulianCenturies(jd)
        return Moon.horizontalParallax(t)
    }

    static distance(jd) {
        const t = toJulianCenturies(jd)
        const theta = Moon.horizontalParallax(t)
        return EARTH_RADIUS_KM / sin(theta)
    }

    // t in julian centuries
    static eclipticLongBase(t) {
        return 218.32 + 481267.881 * t +
             6.29 * sin( 477198.87 * t + 135.0) +
            -1.27 * sin(-413335.36 * t + 259.3) +
             0.66 * sin( 890534.22 * t + 235.7) +
             0.21 * sin( 954397.74 * t + 269.9) +
            -0.19 * sin(  35999.05 * t + 357.5) +
            -0.11 * sin( 966404.03 * t + 186.6)
    }

    static eclipticLatBase(t) {
        const sinConstants = [
            [5.13,  483202.02,  93.3],
            [0.28,  960400.89,  228.2],
            [-0.28, 6003.15,    318.3],
            [-0.17, -407332.21, 217.6]
        ]
        return sum(sinConstants.map(([a, b, c]) => a * sin(b * t + c)))
    }

    static horizontalParallax(t) {
        const cosConstants = [
            [0.0518,    477198.85,  134.9],
            [0.0095,    -413335.38, 259.2],
            [0.0078,    890534.23,  235.7],
            [0.0028,    954397.7,   269.9]
        ] 

        return 0.9508 + sum(cosConstants.map(([a, b, c]) => a * cos(b * t + c)))
    }

    static precessionCorrection(t, lat, long) {
        const a = 1.396971 * t + 0.0003086 * t * t
        const b = 0.013056 * t - 0.0000092 * t * t
        const c = 5.12362 - 1.155358 * t - 0.0001964 * t * t

        const latPre = lat - b * sin(long + c)
        const longPre = long - a + b * cos(long + c) * tan(latPre)
        return [latPre, longPre]
    }
}

function sum(arr) {
    let result = 0;
    for (let a of arr) {
        result += a;
    }
    return result;
}



export class Quaternions {
    // internal [s, v] - external [v, s]
    static toInternalQuat(q) {
        return [q[3], q[0], q[1], q[2]]
    }

    static rotate(vector, quaternion) {
        const quatVector = [0].concat(vector);
        return Quaternions.multiply(
            Quaternions.multiply(quaternion, quatVector),
            Quaternions.inverse(quaternion)
        );
    }

    static squaredNorm(q) {
        return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]
    }

    static multiply(q, r) {
        return [
            r[0] * q[0] - r[1] * q[1] - r[2] * q[2] - r[3] * q[3],
            r[0] * q[1] + r[1] * q[0] - r[2] * q[3] + r[3] * q[2],
            r[0] * q[2] + r[1] * q[3] + r[2] * q[0] - r[3] * q[1],
            r[0] * q[3] - r[1] * q[2] + r[2] * q[1] + r[3] * q[0],
        ];
    }

    static inverse(q) {
        const sn =  Quaternions.squaredNorm(q)
        return [q[0], -q[1], -q[2], -q[3]].map(a => a * 1.0 / sn)
    }
}

// https://www.w3.org/TR/orientation-event/#biblio-eulerangles
function getQuaternion(alpha, beta, gamma) {

  var _x = beta  || 0; // beta value
  var _y = gamma || 0; // gamma value
  var _z = alpha || 0; // alpha value

  var cX = cos( _x/2 );
  var cY = cos( _y/2 );
  var cZ = cos( _z/2 );
  var sX = sin( _x/2 );
  var sY = sin( _y/2 );
  var sZ = sin( _z/2 );

  //
  // ZXY quaternion construction.
  //

  var w = cX * cY * cZ - sX * sY * sZ;
  var x = sX * cY * cZ - cX * sY * sZ;
  var y = cX * sY * cZ + sX * cY * sZ;
  var z = cX * cY * sZ + sX * sY * cZ;

  return [ w, x, y, z ];
}

export function computeAltAzFromQuat(sensorQuaternion) {
    const deviceOriginVector = [0, 0, -1]
    const quaternion = Quaternions.toInternalQuat(sensorQuaternion)
    const directionVec = Quaternions.rotate(deviceOriginVector, quaternion)
    const altitude = toAltitude(directionVec)
    const azimuth = toAzimuth(directionVec)
    return { altitude, azimuth }
}
function toAltitude(vector4) {
    // vector comes from a quaternion ... can throw away scalar 
    const [_, x, y, z] = vector4;
    const vecLenOnXYPlane = Math.sqrt(x**2 + y**2)
    const altitude = atan(z / vecLenOnXYPlane)
    return altitude
}

function toAzimuth(vector4) {
    // vector comes from a quaternion ... can throw away scalar
    const [_, x, y, z] = vector4;

    // [PI, -PI] - positive ccw from east
    const thetaPiMax = -Math.atan2(y, x)

    // [0, 2PI] - positive ccw from east
    const theta2PiMax = thetaPiMax < 0 ? 2 * PI + thetaPiMax : thetaPiMax

    // [0, 2PI] - positive ccw from north
    const thetaFromNorth = (theta2PiMax + PI / 2) % (2 * PI)

    return degree(thetaFromNorth)
}








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
const raToLong = (jd, ra) => (280.1470 + 360.9856235 * (jd - j2000jd) - ra) % 360
const eqCelestialToLatLong = (jd, ra, dec) => ({ lat: dec, long: raToLong(jd, ra) })


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
function getAltAz(jd, userLoc, celestial, print=false) {
    userLoc = toLatLongWest(userLoc) 
    const latLongUnder = eqCelestialToLatLong(jd, celestial.ra, celestial.dec)

    const angleDist = haversineDist(latLongUnder, userLoc)

    if (print) {
        console.log('jd', jd)
        console.log('userLoc', userLoc)
        console.log('celestial', celestial)
        console.log('latLongUnder', latLongUnder)
        console.log('angleDist', angleDist)
    }
    const altitude = 90 - angleDist
    const azimuth = bearing(userLoc, latLongUnder)
    return { altitude, azimuth } 
}

// 0 -> 360 from top to 0->90, 0 to -90
function azToTheta(azimuth) {
    let theta = azimuth
    theta = (theta - 90) % 360
    return theta <= 180 ? -theta : 360 -theta
}

function to3Vec(alt, az, length) {
    const z = length * sin(alt)
    const b = length * cos(alt)
    const theta = azToTheta(az)
    let res = [b*cos(theta), b*sin(theta), z] 
    //console.log(alt, az, length, res)
    return res
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

console.log(getAltAz(toJd(Date.now()), austin, starMirach)) 


function getPosition() {
    return new Promise((resolve, reject) => navigator.geolocation.getCurrentPosition(resolve, reject))
}


const userLoc = await getPosition()
const userLatLong = { lat: userLoc.coords.latitude, long: userLoc.coords.longitude }


const stars = await fetch('./data/stars_vis.json').then(response => response.json())

/*const stars = [
    { name: 'red', mag: 0, dec: 90, ra: 0 }, // np
    { name: 'green', mag: 0, dec: -90, ra: 0 }, //sp
    { name: 'blue', mag: 0, dec: 0, ra: 0 },        //vernal
    { name: 'yellow', mag: 0, dec: 0, ra: 12 }, // op vernal
    { name: 'orange', mag: 0, dec: 0, ra: 6 }, 
    { name: 'purple', mag: 0, dec: 0, ra: 18 } 
]*/

const EARTH_OBLIQUITY = 23.4393 // epsilon


function rightAscension(eclip) {
    const l = cos(eclip.lat) * cos(eclip.long)
    const m = 0.9175 * cos(eclip.lat) * sin(eclip.long) - 0.3978 * sin(eclip.lat)    
    const newRes = atan(m/l) 
    return atan2(sin(eclip.long) * cos(EARTH_OBLIQUITY) - tan(eclip.lat) * sin(EARTH_OBLIQUITY), cos(eclip.long))
    return newRes
}

function declination(eclip) {
    return asin(sin(eclip.lat) * cos(EARTH_OBLIQUITY) + cos(eclip.lat) * sin(EARTH_OBLIQUITY) * sin(eclip.long))
}



// modified to output ra/dec in equitorial rather than 
function eclipticToEquitorial(eclipLatLong, jd) {
    const ra = rightAscension(eclipLatLong)
    const dec = declination(eclipLatLong)
    return { ra, dec } 
}



console.log(stars)

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");


canvas.width  = canvas.clientWidth;
canvas.height = canvas.clientHeight;
const width = canvas.width
const height = canvas.height

console.log('dims', width, height)

// make sphere such that height of phone covers 90 of sphere
const radius = height / Math.sqrt(2) 
console.log('radius', radius)


const starAzFrame3Vec = [
    { name: 'red', mag: 0, x1: 0, y1: 0, z1: radius }, // np
    { name: 'green', mag: 0,  x1: 0, y1: 0, z1: -radius }, //sp
    { name: 'blue', mag: 0, x1:0, y1: radius, z1: 0},        //vernal
    { name: 'yellow', mag: 0, x1:0, y1: -radius, z1: 0}, // op vernal
    { name: 'orange', mag: 0, x1:radius, y1: 0, z1: 0 }, 
    { name: 'purple', mag: 0, x1:-radius, y1: 0, z1:0} 
]

const minMag = 4.5


function drawStar(name, mag, xCanvas, yCanvas) {
    ctx.beginPath();
    // brightest star −1.46 (sirius
    // least bright mag shown 
    //
    let size = -mag + minMag // between 0 and 8.46
    
    const magRange = 1.46 + minMag
    const maxSize = 12
    const minSize = 1
    size = (size / magRange) * (maxSize - minSize) + minSize

    ctx.arc(xCanvas, yCanvas, size, 0, 2 * Math.PI);
    ctx.font = "15pt bold";
    ctx.textAlign = "center";
    ctx.fillText(name, xCanvas, yCanvas + 30);
    ctx.fillStyle = 'white';
    ctx.fill();
}


function drawMoon(xCanvas, yCanvas) {
    ctx.beginPath();
    ctx.arc(xCanvas, yCanvas, 50, 0, 2 * Math.PI);
    ctx.font = "15pt bold";
    ctx.textAlign = "center";
    ctx.fillText("Moon", xCanvas, yCanvas + 30);
    ctx.fillStyle = 'grey';
    ctx.fill();
}

const options = { frequency: 30, referenceFrame: "device" };
const sensor = new AbsoluteOrientationSensor(options);
sensor.start();
sensor.addEventListener("reading", () => {

    //console.log(sensor.quaternion)
    const orientQuat = Quaternions.toInternalQuat(sensor.quaternion)
    const inverseOrientQuat = Quaternions.inverse(orientQuat)

    //console.log(orientQuat )
    const jd = toJd(Date.now())
    //console.log(width, height)

    const totalStars = stars.length
    let renderStars = 0
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    //for (const star of starAzFrame3Vec) {
    //

    const xMax = width / 2, xMin = -width /2, yMax = height / 2, yMin = -height /2
    for (const star of stars) {
        let { name, mag, ra, dec } = star
        // celestial ra from hvg dataset is in hours, so times 15 to get degrees
        ra *= 15

        const { altitude: starAlt, azimuth: starAz} = getAltAz(jd, userLatLong, { ra, dec })
        const star3Vec = to3Vec(starAlt, starAz, radius)
        const [_, x, y, z] = Quaternions.rotate(star3Vec, inverseOrientQuat)
        const inFrame = z < 0 && xMin <= x && x <= xMax && yMin <= y && y <= yMax
        if (mag < minMag && inFrame) {
            renderStars += 1
            const xCanvas = x + width / 2
            const yCanvas = -y + height / 2
            drawStar(name, mag, xCanvas, yCanvas)
        }
    }

    console.log('putz')

    const moonEclipLL = Moon.eclipLatLong(jd)
    const moonLoc = eclipticToEquitorial(moonEclipLL, jd)
    const { altitude: moonAlt, azimuth: moonAz} = getAltAz(jd, userLatLong, moonLoc, true)
    const moon3Vec = to3Vec(moonAlt, moonAz, radius)
    const [_, mx, my, mz] = Quaternions.rotate(moon3Vec, inverseOrientQuat)


    console.log('moonEclipLL', moonEclipLL)
    console.log('moonLoc', moonLoc)
    console.log('moon alt/az', moonAlt, moonAz)
    console.log('moon3Vec', moon3Vec)

    if (mz < 0) {
        const xmCanvas = mx + width / 2
        const ymCanvas = -my + height / 2
    console.log('moon canvas', xmCanvas, ymCanvas)
        drawMoon(xmCanvas, ymCanvas)
    }




    console.log('totalStars', totalStars)
    console.log('renderStars', renderStars)

    //console.log('phone alt/az', altAz)
});
sensor.addEventListener("error", (error) => {
    console.log(error)
  if (event.error.name === "NotReadableError") {
    console.log("Sensor is not available.");
  }
});



