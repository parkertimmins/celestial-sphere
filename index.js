//----------------------- General Math and utilities -----------------------
const mod = (m, n) => ((m%n)+n)%n 
const rad = (deg) => deg * Math.PI / 180
const degree = (radian) => radian * (180 / Math.PI)
const sin = (deg) => Math.sin(rad(deg))
const cos = (deg) => Math.cos(rad(deg))
const tan = (deg) => Math.tan(rad(deg))
const acos = (x) => degree(Math.acos(x))
const asin = (x) => degree(Math.asin(x))
const atan2 = (x, y) => degree(Math.atan2(x, y))
const sum = (arr) => arr.reduce((a, b) => a+b, 0)
const squaredNorm = (v) => sum(v.map(e => e*e))
const scalarMult = (a, v) => v.map(x => a*x)
const euclideanDist = (x1, y1, x2, y2) => Math.sqrt((x2-x1)**2 + (y2-y1)**2)

class Quaternions {
    static fromAngleAxis(angle, axis) {
        return [cos(angle/2)].concat(scalarMult(sin(angle/2), axis))
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
        return scalarMult(1/squaredNorm(q), [q[0], -q[1], -q[2], -q[3]])
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

// https://en.wikipedia.org/wiki/Exponential_smoothing
const expAvgFilter = (a = 0.5) => {
    return {
        x: null,
        update(xt) { 
            this.x = (this.x === null ? xt : (a * xt + (1 - a) * this.x)) 
        },
        get value () { return this.x }
    }
}


//----------------------- Julian date -----------------------
// https://www.aa.quae.nl/en/reken/zonpositie.html
const millisPerDay = 1000 * 60 * 60 * 24
const daysPerCentury = 36525
const j2000Date = Date.UTC(2000, 0, 1, 12, 0, 0)
const j2000jd = 2451545
const toJulianCenturies = (jd) =>  (jd - j2000jd) / daysPerCentury
const toJd = (date) => (date - j2000Date) / millisPerDay + j2000jd


//----------------------- Compute ecliptic coordinates of objects that move against celestial sphere -----------------------
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
    const [x, y, z] = [xp - xe, yp - ye, zp - ze]
    const delta = Math.sqrt(x*x + y*y + z*z)
    const lambda = atan2(y, x)
    const beta = asin(z/delta)
    return { long: lambda, lat: beta } 
}

//----------------------- Conversion from Ecliptic to Equitorial -----------------------

// change to coordinates on a sphere shifted by angle along with 90 degree meridian
function sphereCoordTransform({lat, long}, angle) {
    const longAngle = atan2(sin(long) * cos(angle) - tan(lat) * sin(angle), cos(long))
    const latAngle = asin(sin(lat) * cos(angle) + cos(lat) * sin(angle) * sin(long))
    return {lat: latAngle, long: longAngle }
}

const EARTH_OBLIQUITY = 23.4393
function eclipticToEquitorial(eclip) {
    const {lat, long} = sphereCoordTransform(eclip, EARTH_OBLIQUITY)
    return { ra: long, dec: lat }
}

//----------------------- Conversion from Equitorial to Altitude/Azimuth (horizontal coorinates) -----------------------
const raToLong = (jd, ra) => (280.1470 + 360.9856235 * (jd - j2000jd) - ra) % 360
const equitorialToLatLong = (jd, {ra, dec}) => ({ lat: dec, long: raToLong(jd, ra) })

// computer altitude/azimuth from ra/dec of a celestial object at current time, and user location in lat long
// celestial object is in equitorial frame
// assume parallax is 0, eg at inifinite distance
function getAltAz(user, celestial) {
    const angle = 90 - user.lat
    const point = { lat: celestial.lat, long: celestial.long - user.long + 90 }
    const {lat, long} = sphereCoordTransform(point, angle)
    return { altitude: lat, azimuth: mod(long + 90, 360) } 
}

// theta: west=0, positive to north, negative to south, azimuth: north=0, increase ccw, all positive
const thetaToAz = (theta) => mod(-theta + 90, 360)

function azToTheta(azimuth) {
  const theta = (azimuth - 90) % 360
  return theta <= 180 ? -theta : 360 - theta
}

// long: negative to west, positive to east, longWest: increase to west, all positive
const toLongWest = long => mod(-long, 360)


// to vectors on unit sphere
function to3Vec(alt, az) {
    const b = cos(alt)
    const theta = azToTheta(az)
    return [b*cos(theta), b*sin(theta), sin(alt)] 
}

//----------------------- Rendering/UI code below -----------------------

function getPosition() {
    return new Promise((resolve, reject) => navigator.geolocation.getCurrentPosition(resolve, reject))
}
const loadImage = (url) => new Promise((resolve, reject) => {
  const img = new Image();
  img.addEventListener('load', () => resolve(img));
  img.src = url;
});


const toPixelSize = (deg, st) => st.bounds.sphereToPixScale * rad(deg)
const drawImgCentered = (ctx, img, x, y, size) => ctx.drawImage(img, x-size/2, y-size/2, size, size)

function addTitle(st, ctx, x, y, text, color, pixSize) {
    ctx.font = '22pt Calibri';
    ctx.textAlign = "center";
    ctx.fillStyle = color;
    const textPixOffset = pixSize / 2 + toPixelSize(1, st);
    ctx.fillText(text, x, y + textPixOffset);
}

function drawStar(st, name, mag, x, y) {
    const percMagRange = (-mag + minVisibleMag) / magRange // flip [-1.46, 4.5] and map to [0, 1]
    const imgRadius = percMagRange * (maxStarSize(st) - minStarSize(st)) + minStarSize(st)
    ctx.beginPath();
    ctx.arc(x, y, imgRadius, 0, 2 * Math.PI);
    ctx.fillStyle = 'white';
    ctx.fill();
    addTitle(st, ctx, x, y, name, 'white', imgRadius);
}

function drawMoon(st, x, y) {
    const pixSize = toPixelSize(2.5, st)
    drawImgCentered(ctx, images.Moon, x, y, pixSize)
    addTitle(st, ctx, x, y, 'Moon', 'lightgreen', pixSize);
}

function drawSun(st, x, y) {
    const pixSize = toPixelSize(2.5, st)
    drawImgCentered(ctx, images.Sun, x, y, pixSize) 
    addTitle(st, ctx, x, y, 'Sun', 'white', pixSize);
}

function drawPlanet(st, p, x, y) {
    const pixSize = toPixelSize(p.imgSize, st)
    drawImgCentered(ctx, images[p.name], x, y, pixSize) 
    addTitle(st, ctx, x, y, p.name, 'lightgreen', pixSize);
}

function computeBounds(longVisAngle) {
    const distToPlane = cos(longVisAngle / 2)
    const hSphere = 2 * sin(longVisAngle / 2)
    const wSphere = canvas.width * (hSphere / canvas.height)
    const sphereToPixScale = canvas.height / hSphere
    const xBase = -wSphere/2, yBase = -hSphere/2
    // more min/max larger than screen so off screen rendered cleanly as enters frame
    const xMax = wSphere, xMin = -wSphere, yMax = hSphere, yMin = -hSphere
    return { xBase, yBase, xMax, xMin, yMax, yMin, sphereToPixScale, distToPlane }
}

function toCanvasCoords(jd, ra, dec, inverseOrientQuat, bounds, userLatLong) {
    const celestialLatLong = equitorialToLatLong(jd, {ra, dec})
    const { altitude, azimuth } = getAltAz(userLatLong, celestialLatLong)
    const vecOnSphere = to3Vec(altitude, azimuth)
    const rotVecOnSphere = Quaternions.rotate(vecOnSphere, inverseOrientQuat).slice(1)
    if (rotVecOnSphere[2] > 0) {
        return null
    }

    //https://math.stackexchange.com/questions/3412199/how-to-calculate-the-intersection-point-of-a-vector-and-a-plane-defined-as-a-poi
    const [x, y, z] = scalarMult(-bounds.distToPlane / rotVecOnSphere[2], rotVecOnSphere)
    if (bounds.xMin <= x && x <= bounds.xMax && bounds.yMin <= y && y <= bounds.yMax) {
        const xPixOff = (x - bounds.xBase) * bounds.sphereToPixScale 
        const yPixOff = (y - bounds.yBase) * bounds.sphereToPixScale 
        return [xPixOff, canvas.height - yPixOff]
    } 
    return null;
}

function addObject(st, inverseOrientQuat, jd, eclipLatLong, drawFunc) {
    const {ra, dec} = eclipticToEquitorial(eclipLatLong)
    const coords = toCanvasCoords(jd, ra, dec, inverseOrientQuat, st.bounds, st.userLatLong)
    if (coords !== null) {
        drawFunc(st, ...coords)
    }
}

function render(st, ctx, canvas) {
    const inverseOrientQuat = Quaternions.inverse(st.orientQuat)
    const jd = toJd(Date.now())
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    for (const star of stars) {
        let { name, mag, ra, dec } = star
        if (mag < minVisibleMag) {
            const coords = toCanvasCoords(jd, ra, dec, inverseOrientQuat, st.bounds, st.userLatLong)
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

const brighestStarMag = -1.46           // sirius
const minVisibleMag = 6               // dimmest magnitude shown
const magRange = -brighestStarMag + minVisibleMag

// Start loading files for stars, planet data, and images
const promises = []
promises.push(fetch('./data/stars/yale_stars.json').then(response => response.json()))
promises.push(fetch('./data/planets.json').then(response => response.json()))
const objects = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
objects.forEach(obj => promises.push(loadImage(`./images/icons/${obj}.png`)))


// Wait till files load 
const resolved = await Promise.all(promises)
const stars = resolved[0].map(s => ({ra:s[0], dec:s[1], mag:s[2], name:(s.length===4 ? s[3] : '')}))
const planets = resolved[1]
const earth = planets.filter((p) => p.name === "Earth")[0]
const images = Object.fromEntries(objects.map((name, i) => [name, resolved[i+2]]))


// Global mutable state
const state = {} 
// degrees of sky covered by long edge of screen 
state.longVisAngle = 90
// quaternion describing current phone orientation
state.orientQuat = [0, 0, 0, 0]
// iPhone saves difference between relative north and absolute north, filtered to avoid jumps
state.bearingDiffFilter = expAvgFilter(0.05)
// data derived from longVisAngle, held in state for efficiency
state.bounds = computeBounds(state.longVisAngle)
// cache of pointer events for zooming 
state.evCache = {}
// lat/long of user
state.userLatLong = {lat:0, long: 0}


// All functions that access state object directly are below this point
function iosRenderOnOrientChange() {
    if (typeof DeviceOrientationEvent.requestPermission === 'function') {
      DeviceOrientationEvent.requestPermission(true)
        .then(permissionState => {
          if (permissionState === 'granted') {

            window.addEventListener('deviceorientation', () => {
                if (event.webkitCompassHeading === 0 && state.bearingDiffFilter.value === null) {
                    console.log('skipping initial orientation');
                    return;
                }

                const relativeQuat = Quaternions.fromAngles(event.alpha, event.beta, event.gamma)
                const phoneNorth = [0, 1, 0] 
                const phoneBack = [0, 0, -1] 
                const northRotated = Quaternions.rotate(phoneNorth, relativeQuat).slice(1)
                const backRotated = Quaternions.rotate(phoneBack, relativeQuat).slice(1)
                const lenProjOnXy = Math.sqrt(northRotated[0]**2 + northRotated[1]**2)
                const northRotatedAngleFromHorizon = Math.abs(atan2(northRotated[2], lenProjOnXy))

                if (backRotated[2] < 0 && northRotatedAngleFromHorizon < 80) {
                    const thetaRelativeNorth = atan2(northRotated[1], northRotated[0])
                    const bearingRelativeNorth = thetaToAz(thetaRelativeNorth)
                    const bearingDiff = mod(event.webkitCompassHeading - bearingRelativeNorth, 360)
                    state.bearingDiffFilter.update(bearingDiff)
                }
               
                const noHeadingAlert = document.getElementById("alert-container")
                noHeadingAlert.style.display = state.bearingDiffFilter.value === null ? 'flex' : 'none';
                const angle = state.bearingDiffFilter.value || 0;
                const northOffsetQuat = Quaternions.fromAngleAxis(angle, [0, 0, -1])
                state.orientQuat = Quaternions.multiply(northOffsetQuat, relativeQuat)
                
                render(state, ctx, canvas)
            }, true);
          }
        })
        .catch(console.error);
    } 
}

function androidRenderOnOrientChange() {
    window.addEventListener("deviceorientationabsolute", (event) => {
        console.log(event);
        state.orientQuat = Quaternions.fromAngles(event.alpha, event.beta, event.gamma)
        render(state, ctx, canvas)
    });
}

function setCanvasSize() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
    state.bounds = computeBounds(state.longVisAngle)
}

const isIOS = () => /(iPad|iPhone)/g.test(navigator.userAgent)
const allowButton = document.getElementById("request-perms")
allowButton.onclick = () => {
    document.getElementById("request-perms-pane").style.display = 'none';
    document.getElementById("star-pane").style.display = 'block';
    setCanvasSize();

    // start gps 
    navigator.geolocation.watchPosition(pos => {
        state.userLatLong = { lat: pos.coords.latitude, long: toLongWest(pos.coords.longitude) }
    }, console.log);
  
    if (isIOS()) {
        iosRenderOnOrientChange() 
    } else {
        androidRenderOnOrientChange() 
    }   
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
        state.longVisAngle = Math.min(90, prevDegDist * (canvas.height / currPixDist)) 
        state.bounds = computeBounds(state.longVisAngle)
        render(state, ctx, canvas)
    } else {
        state.evCache[ev.pointerId] = ev
    }
}

