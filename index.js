
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
      atan2 = (x, y) => degree(Math.atan2(x, y));

// Simple vector functions
const sum = (arr) => arr.reduce((a, b) => a+b, 0)
const squaredNorm = (v) => sum(v.map(e => e*e))
const norm = (v) => Math.sqrt(squaredNorm(v))
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
const dot = (a, b) => sum(a.map((x, i) => x*b[i]))
const scalarMult = (a, v) => v.map((x, i) => a * x)
const euclideanDist = (x1, y1, x2, y2) => Math.sqrt((x2-x1)**2 + (y2-y1)**2)

// Originally from most recent Astronomical Almanac, at least 1999-2015
// An Alternative Lunar Ephemeris Model
// https://caps.gsfc.nasa.gov/simpson/pubs/slunar.pdf
// Ignore parallax
function moonEclipLatLong(jd) {
    const t = toJulianCenturies(jd)
    const long = 218.32 + 481267.881 * t +
         6.29 * sin( 477198.87 * t + 135.0) +
        -1.27 * sin(-413335.36 * t + 259.3) +
         0.66 * sin( 890534.22 * t + 235.7) +
         0.21 * sin( 954397.74 * t + 269.9) +
        -0.19 * sin(  35999.05 * t + 357.5) +
        -0.11 * sin( 966404.03 * t + 186.6)
    
    const lat = 5.13 * sin( 483202.02 * t + 93.3) +
        0.28 * sin( 960400.89 * t +  228.2) +
       -0.28 * sin( 6003.15 * t +    318.3) +
       -0.17 * sin( -407332.21 * t + 217.6)
    
    return { lat, long }
}

class Quaternions {
    static fromAngleAxis(angle, axis3Vec) {
        const sinAngle2 = sin(angle/2)
        return [cos(angle/2)].concat(axis3Vec.map(a => a*sinAngle2))
    }

    // internal [s, v] - external [v, s]
    static toInternalQuat(q) {
        return [q[3], q[0], q[1], q[2]]
    }

    static rotate(vector, quaternion) {
        return Quaternions.multiply(
            Quaternions.multiply(quaternion, [0].concat(vector)),
            Quaternions.inverse(quaternion)
        );
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
        const sn =  squaredNorm(q)
        return [q[0], -q[1], -q[2], -q[3]].map(a => a / sn)
    }
    
    // https://w3c.github.io/deviceorientation/#deviceorientation (Example 11)
    static fromAngles(alpha, beta, gamma) {
      var _x = beta  || 0; 
      var _y = gamma || 0; 
      var _z = alpha || 0;

      var cX = cos( _x/2 );
      var cY = cos( _y/2 );
      var cZ = cos( _z/2 );
      var sX = sin( _x/2 );
      var sY = sin( _y/2 );
      var sZ = sin( _z/2 );

      // ZXY quaternion construction.
      var w = cX * cY * cZ - sX * sY * sZ;
      var x = sX * cY * cZ - cX * sY * sZ;
      var y = cX * sY * cZ + sX * cY * sZ;
      var z = cX * cY * sZ + sX * sY * cZ;

      return [ w, x, y, z ];
    }
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
const raToLong = (jd, ra) => (280.1470 + 360.9856235 * (jd - j2000jd) - ra) % 360

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
    const latLongUnder = { lat: celestial.dec, long: raToLong(jd, celestial.ra) }
    const angleDist = haversineDist(latLongUnder, userLoc)
    const altitude = 90 - angleDist
    const azimuth = bearing(userLoc, latLongUnder)
    return { altitude, azimuth } 
}

const thetaToAz = (theta) => mod(-theta + 90, 360)

// 0 -> 360 from top to 0->90, 0 to -90
function azToTheta(azimuth) {
    let theta = (azimuth - 90) % 360
    return theta <= 180 ? -theta : 360 - theta
}

function toLatLongWest(latLong) {
    const {lat, long} = latLong
    const longWest = long < 0 ? -long : 360 - long
    return {lat: lat, long: longWest } 
}

// to vectors on unit sphere
function to3Vec(alt, az) {
    const b = cos(alt)
    const theta = azToTheta(az)
    return [b*cos(theta), b*sin(theta), sin(alt)] 
}

function getPosition() {
    return new Promise((resolve, reject) => navigator.geolocation.getCurrentPosition(resolve, reject))
}
const loadImage = (url) => new Promise((resolve, reject) => {
  const img = new Image();
  img.addEventListener('load', () => resolve(img));
  img.src = url;
});


const EARTH_OBLIQUITY = 23.4393
const rightAscension = (eclip) => atan2(sin(eclip.long) * cos(EARTH_OBLIQUITY) - tan(eclip.lat) * sin(EARTH_OBLIQUITY), cos(eclip.long))
const declination = (eclip) => asin(sin(eclip.lat) * cos(EARTH_OBLIQUITY) + cos(eclip.lat) * sin(EARTH_OBLIQUITY) * sin(eclip.long))
const eclipticToEquitorial = (eclip) => ({ ra: rightAscension(eclip), dec: declination(eclip) })

//https://www.aa.quae.nl/en/reken/hemelpositie.html
//https://www.aa.quae.nl/en/reken/zonpositie.html for kepler approximation
function rectHelioEcliptical(jd, p) {
    const d = jd
    const d0 = j2000jd
    const M = (p.M0 + p.n * (d - d0)) % 360
    const C = p.C1 * sin(M) + p.C2 * sin(2 * M) + p.C3 * sin(3 * M) * p.C4 * sin(4 * M) + p.C5 * sin(5 * M) + p.C6 * sin(6 * M)
    const nu = (M + C) % 360
    const r = p.a * (1 - p.e*p.e) / (1 + p.e * cos(nu))
    const x = r * (cos(p.OMEGA) * cos(p.omega + nu) - sin(p.OMEGA) * cos(p.i) * sin(p.omega + nu))
    const y = r * (sin(p.OMEGA) * cos(p.omega + nu) + cos(p.OMEGA) * cos(p.i) * sin(p.omega + nu))
    const z = r * sin(p.i) * sin(p.omega + nu)
    return [x, y, z]    
}

function sunEclipLatLong(jd, earth) {
    const [x, y, z] = rectHelioEcliptical(jd, earth) 
    // add 180 since sun from earth rather than earth from sun 
    return { long: mod(atan2(y, x) + 180, 360), lat: 0 } 
}

// https://www.aa.quae.nl/en/reken/zonpositie.html
function planetEclipLatLong(jd, p, earth) {
    const [xp, yp, zp] = rectHelioEcliptical(jd, p) 
    const [xe, ye, ze] = rectHelioEcliptical(jd, earth) 
    const x = xp - xe
    const y = yp - ye
    const z = zp - ze
    const delta = Math.sqrt(x*x + y*y + z*z)
    const lambda = atan2(y, x)
    const beta = asin(z/delta)
    return { long: lambda, lat: beta } 
}

function toPixelSize(deg, st) {
    return st.bounds.sphereToPixScale * rad(deg)
}

const drawImgCentered = (ctx, img, x, y, size) => ctx.drawImage(img, x-size/2, y-size/2, size, size)

function addTitle(st, ctx, x, y, text, color, font, pixSize) {
    ctx.font = font;
    ctx.textAlign = "center";
    ctx.fillStyle = color;
    const textPixOffset = pixSize / 2 + toPixelSize(1, st)
    ctx.fillText(text, x, y + textPixOffset);
}

function drawStar(st, name, mag, x, y) {
    const percMagRange = (-mag + minVisibleMag) / magRange // flip [-1.46, 4.5] and map to [0, 1]
    const imgRadius = percMagRange * (maxStarSize(st) - minStarSize(st)) + minStarSize(st)
    ctx.beginPath();
    ctx.arc(x, y, imgRadius, 0, 2 * Math.PI);
    ctx.fillStyle = 'white';
    ctx.fill();
    addTitle(st, ctx, x, y, name, 'white', "20pt bold", imgRadius);
}

function drawMoon(st, x, y) {
    const pixSize = toPixelSize(2.5, st)
    drawImgCentered(ctx, images.Moon, x, y, pixSize)
    addTitle(st, ctx, x, y, 'Moon', 'lightgreen', "20pt bold", pixSize);
}

function drawSun(st, x, y) {
    const pixSize = toPixelSize(2.5, st)
    drawImgCentered(ctx, images.Sun, x, y, pixSize) 
    addTitle(st, ctx, x, y, 'Sun', 'white', "20pt bold", pixSize);
}

function drawPlanet(st, p, x, y) {
    const pixSize = toPixelSize(p.imgSize, st)
    drawImgCentered(ctx, images[p.name], x, y, pixSize) 
    addTitle(st, ctx, x, y, p.name, 'lightgreen', "20pt bold", pixSize);
}

function computeBounds(longVisAngle) {
    const distToPlane = cos(longVisAngle / 2)
    const hSphere = 2 * sin(longVisAngle / 2)
    const heightToWidthRatio = canvas.height / canvas.width
    const wSphere = hSphere / heightToWidthRatio
    const sphereToPixScale = canvas.height / hSphere
    const xBase = -wSphere/2, yBase = -hSphere/2
    // more min/max larger than screen so off screen rendered cleanly as enters frame
    const xMax = wSphere, xMin = -wSphere, yMax = hSphere, yMin = -hSphere
    return { xBase, yBase, xMax, xMin, yMax, yMin, sphereToPixScale, distToPlane }
}

function toCanvasCoords(jd, ra, dec, inverseOrientQuat, bounds) {
    const { altitude, azimuth} = getAltAz(jd, userLatLong, { ra, dec })
    const vecOnSphere = to3Vec(altitude, azimuth)
    const rotVecOnSphere = Quaternions.rotate(vecOnSphere, inverseOrientQuat).slice(1)
    if (rotVecOnSphere[2] > 0) {
        return null
    }

    //https://math.stackexchange.com/questions/3412199/how-to-calculate-the-intersection-point-of-a-vector-and-a-plane-defined-as-a-poi
    const [x, y, z] = scalarMult(-bounds.distToPlane / rotVecOnSphere[2], rotVecOnSphere)
    const inFrame = bounds.xMin <= x && x <= bounds.xMax && bounds.yMin <= y && y <= bounds.yMax
    if (inFrame) {
        const xPixOff = (x - bounds.xBase) * bounds.sphereToPixScale 
        const yPixOff = (y - bounds.yBase) * bounds.sphereToPixScale 
        return [xPixOff, canvas.height - yPixOff]
    } else {
        return null;
    }
}

function addObject(st, inverseOrientQuat, jd, eclipLatLong, drawFunc) {
    const {ra, dec} = eclipticToEquitorial(eclipLatLong)
    const coords = toCanvasCoords(jd, ra, dec, inverseOrientQuat, st.bounds)
    if (coords !== null) {
        drawFunc(st, ...coords)
    }
}

function render(st) {
    const inverseOrientQuat = Quaternions.inverse(st.orientQuat)
    const jd = toJd(Date.now())
    ctx.clearRect(0, 0, canvas.width, canvas.height); // TODO canvas ctx

    for (const star of stars) {
        let { name, mag, ra, dec } = star
        if (mag < minVisibleMag) {
            const coords = toCanvasCoords(jd, ra, dec, inverseOrientQuat, st.bounds)
            if (coords !== null) {
                drawStar(st, name, mag, coords[0], coords[1])
            }
        }
    }
    
    addObject(st, inverseOrientQuat, jd, moonEclipLatLong(jd), drawMoon);
    addObject(st, inverseOrientQuat, jd, sunEclipLatLong(jd, earth), drawSun);
    for (const p of planets) {
        if (p.name !== "Earth") {
            addObject(st, inverseOrientQuat, jd, planetEclipLatLong(jd, p, earth), (st, x, y) => drawPlanet(st, p, x, y));
        }
    }
}


const maxStarSize = (st) => toPixelSize(0.3, st) 
const minStarSize = (st) => toPixelSize(0.01, st)

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
canvas.width = window.innerWidth;
canvas.height = window.innerHeight;

const brighestStarMag = -1.46           // sirius
const minVisibleMag = 6               // dimmest magnitude shown
const magRange = -brighestStarMag + minVisibleMag

// Start loading files for stars, planet data, and images
const promises = []
promises.push(fetch('./data/stars.json').then(response => response.json()))
promises.push(fetch('./data/planets.json').then(response => response.json()))
const objects = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
objects.forEach(obj => promises.push(loadImage(`./images/icons/${obj}.jpg`)))


// Wait till files load 
const resolved = await Promise.all(promises)
const stars = resolved[0]
const planets = resolved[1]
const earth = planets.filter((p) => p.name === "Earth")[0]
const images = {}
for (let i = 0; i < objects.length; i++) {
    images[objects[i]] = resolved[i+2]
}

const userLoc = await getPosition()
const userLatLong = { lat: userLoc.coords.latitude, long: userLoc.coords.longitude }


// Global mutable state
const state = {} 
// degrees of sky covered by long edge of screen 
state.longVisAngle = 90
// quaternion describing current phone orientation
state.orientQuat = [0, 0, 0, 0]
// iPhone saves offset quaternion between relative north and absolute north
state.northOffsetQuat = null
// data derived from longVisAngle, held in state for efficiency
state.bounds = computeBounds(state.longVisAngle)
// cache of pointer events for zooming 
state.evCache = {}


// All functions that access state object directly are below this point
function iosRenderOnOrientChange() {
    document.getElementById("request-perms").style.display = 'none';
    if (typeof DeviceOrientationEvent.requestPermission === 'function') {
      DeviceOrientationEvent.requestPermission()
        .then(permissionState => {
          if (permissionState === 'granted') {
            window.addEventListener('deviceorientation', () => {
                const relativeQuat = Quaternions.fromAngles(event.alpha, event.beta, event.gamma)
                if (state.northOffsetQuat === null) {
                    const phoneNorth = [0, 1, 0]
                    const northRotated = Quaternions.rotate(phoneNorth, relativeQuat).slice(1)
                    const thetaRelativeNorth = atan2(northRotated[1], northRotated[0])
                    const bearingRelativeNorth = thetaToAz(thetaRelativeNorth)
                    const bearingDiff = event.webkitCompassHeading - bearingRelativeNorth 
                    state.northOffsetQuat = Quaternions.fromAngleAxis(bearingDiff, [0, 0, -1]) 
                }
                state.orientQuat = Quaternions.multiply(state.northOffsetQuat, relativeQuat)
                render(state)
            });
          }
        })
        .catch(console.error);
    } 
}


function androidRenderOnOrientChange() {
    const options = { frequency: 30, referenceFrame: "device" };
    const sensor = new AbsoluteOrientationSensor(options);
    sensor.start();
    sensor.addEventListener("reading", () => {
        state.orientQuat = Quaternions.toInternalQuat(sensor.quaternion)
        render(state)
     });
    sensor.addEventListener("error", (error) => console.log(error));
}

// Render on orientation change
if (/(iPad|iPhone)/g.test(navigator.userAgent)) {
    const allowButton = document.getElementById("request-perms")
    allowButton.style.display = 'block';
    allowButton.onclick = iosRenderOnOrientChange;
} else {
    androidRenderOnOrientChange() 
}

// Render on zoom in/out
// https://developer.mozilla.org/en-US/docs/Web/API/Pointer_events/Pinch_zoom_gestures
canvas.onpointerdown = (ev) => { state.evCache[ev.pointerId] = ev }
canvas.onpointerup = canvas.onpointercancel = canvas.onpointerout = canvas.onpointerleave = (ev) => { delete state.evCache[ev.pointerId] };
canvas.onpointermove = canvas.onpointermove = (ev) => {
    if (Object.keys(state.evCache).length == 2) {
        const [a1, b1] = Object.values(state.evCache)
        const prevPixDist = euclideanDist(a1.clientX, a1.clientY, b1.clientX, b1.clientY)
        const prevDegDist = state.longVisAngle * (prevPixDist / canvas.height)
        state.evCache[ev.pointerId] = ev
        const [a2, b2] = Object.values(state.evCache)
        const currPixDist = euclideanDist(a2.clientX, a2.clientY, b2.clientX, b2.clientY)
        const newLongVisAngle = Math.min(90, prevDegDist * (canvas.height / currPixDist)) 
        state.longVisAngle = newLongVisAngle // TODO global
        state.bounds = computeBounds(state.longVisAngle)
        render(state)
    } else {
        state.evCache[ev.pointerId] = ev
    }
}

