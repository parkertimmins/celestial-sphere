
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
class Moon {
    static eclipLatLong(jd) {
        const t = toJulianCenturies(jd)
        return { long: Moon.eclipticLong(t), lat: Moon.eclipticLat(t) }
    }

    // t in julian centuries
    static eclipticLong(t) {
        return 218.32 + 481267.881 * t +
             6.29 * sin( 477198.87 * t + 135.0) +
            -1.27 * sin(-413335.36 * t + 259.3) +
             0.66 * sin( 890534.22 * t + 235.7) +
             0.21 * sin( 954397.74 * t + 269.9) +
            -0.19 * sin(  35999.05 * t + 357.5) +
            -0.11 * sin( 966404.03 * t + 186.6)
    }

    static eclipticLat(t) {
        const sinConstants = [
            [5.13,  483202.02,  93.3],
            [0.28,  960400.89,  228.2],
            [-0.28, 6003.15,    318.3],
            [-0.17, -407332.21, 217.6]
        ]
        return sum(sinConstants.map(([a, b, c]) => a * sin(b * t + c)))
    }
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
    let theta = azimuth
    theta = (theta - 90) % 360
    return theta <= 180 ? -theta : 360 -theta
}

// to vectors on unit sphere
function to3Vec(alt, az) {
    const b = cos(alt)
    const theta = azToTheta(az)
    return [b*cos(theta), b*sin(theta), sin(alt)] 
}

function toLatLongWest(latLong) {
    const {lat, long} = latLong
    const longWest = long < 0 ? -long : 360 - long
    return {lat: lat, long: longWest } 
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

const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
canvas.width  = canvas.clientWidth;
canvas.height = canvas.clientHeight;

const width = canvas.width
const height = canvas.height
const heightToWidthRatio = height / width




let orientQuat = [0, 1, 1, 1]
let longVisAngle = 90 // start with 1/4 screen visible



function toPixelSize(deg, visAngle) {
    const hSphere = 2 * sin(visAngle / 2)
    const wSphere = hSphere / heightToWidthRatio
    const sphereToPixScale = height / hSphere
    return sphereToPixScale * rad(deg)
}
const maxStarSize = (visAngle) => toPixelSize(0.3, visAngle) 
const minStarSize = (visAngle) => toPixelSize(0.01, visAngle)

const drawImgCentered = (ctx, img, x, y, size) => ctx.drawImage(img, x-size/2, y-size/2, size, size)
const brighestStarMag = -1.46           // sirius
const minVisibleMag = 6               // dimmest magnitude shown
const magRange = -brighestStarMag + minVisibleMag


function addTitle(ctx, x, y, text, color, font, pixSize) {
    ctx.font = font;
    ctx.textAlign = "center";
    ctx.fillStyle = color;
    const textPixOffset = pixSize / 2 + toPixelSize(1, longVisAngle)
    ctx.fillText(text, x, y + textPixOffset);
}

function drawStar(name, mag, x, y) {
    const percMagRange = (-mag + minVisibleMag) / magRange // flip [-1.46, 4.5] and map to [0, 1]
    const imgRadius = percMagRange * (maxStarSize(longVisAngle) - minStarSize(longVisAngle)) + minStarSize(longVisAngle)
    
    ctx.beginPath();
    ctx.arc(x, y, imgRadius, 0, 2 * Math.PI);
    ctx.fillStyle = 'white';
    ctx.fill();
    
    addTitle(ctx, x, y, name, 'white', "20pt bold", imgRadius);
}


function drawMoon(x, y) {
    const pixSize = toPixelSize(2.5, longVisAngle)
    drawImgCentered(ctx, images.Moon, x, y, pixSize)
    addTitle(ctx, x, y, 'Moon', 'lightgreen', "20pt bold", pixSize);
}

function drawSun(x, y) {
    const pixSize = toPixelSize(2.5, longVisAngle)
    drawImgCentered(ctx, images.Sun, x, y, pixSize) 
    addTitle(ctx, x, y, 'Sun', 'white', "20pt bold", pixSize);
}

function drawPlanet(p, x, y) {
    const pixSize = toPixelSize(p.imgSize, longVisAngle)
    drawImgCentered(ctx, images[p.name], x, y, pixSize) 
    addTitle(ctx, x, y, p.name, 'lightgreen', "20pt bold", pixSize);
}

function toCanvasCoords(jd, ra, dec, inverseOrientQuat) {
    const { altitude, azimuth} = getAltAz(jd, userLatLong, { ra, dec })
    const vecOnSphere = to3Vec(altitude, azimuth)
    const rotVecOnSphere = Quaternions.rotate(vecOnSphere, inverseOrientQuat).slice(1)
    if (rotVecOnSphere[2] > 0) {
        return null
    }

    //https://math.stackexchange.com/questions/3412199/how-to-calculate-the-intersection-point-of-a-vector-and-a-plane-defined-as-a-poi
    const distToPlane = cos(longVisAngle / 2)
    const [x, y, z] = scalarMult(-distToPlane / rotVecOnSphere[2], rotVecOnSphere)

    // TODO only recompute on resize
    const hSphere = 2 * sin(longVisAngle / 2)
    const wSphere = hSphere / heightToWidthRatio
    const sphereToPixScale = height / hSphere
    const xMax = wSphere/2, xMin = -wSphere/2, yMax = hSphere/2, yMin = -hSphere/2
    
    const inFrame = z < 0 && xMin <= x && x <= xMax && yMin <= y && y <= yMax
    if (inFrame) {
        const xPixOff = (x - xMin) * sphereToPixScale 
        const yPixOff = (y - yMin) * sphereToPixScale 
        return [xPixOff, height - yPixOff]
    } else {
        return null;
    }
}

// https://developer.mozilla.org/en-US/docs/Web/API/Pointer_events/Pinch_zoom_gestures
class PinchZoom {
    constructor(element) {
        this.evCache = []
    }

    removeEvent(ev) {
        const index = this.evCache.findIndex((cachedEv) => cachedEv.pointerId === ev.pointerId)
        if (index >= 0) {
            this.evCache.splice(index, 1)
        }
    }

    replaceEvent(ev) {
        const index = this.evCache.findIndex((cachedEv) => cachedEv.pointerId === ev.pointerId);
        if (index >= 0) {
            this.evCache[index] = ev;
        }
    }

    onDownHandler(ev) {
        this.evCache.push(ev);
    }

    onUpHandler(ev) {
        this.removeEvent(ev)
    }

    onMoveHandler(ev) {
        if (this.evCache.length == 2) {
            const prevPixDist = euclideanDist(this.evCache[0].clientX, this.evCache[0].clientY, this.evCache[1].clientX, this.evCache[1].clientY)
            const prevDegDist = longVisAngle * (prevPixDist / height)
            this.replaceEvent(ev)
            const currPixDist = euclideanDist(this.evCache[0].clientX, this.evCache[0].clientY, this.evCache[1].clientX, this.evCache[1].clientY)
            const newLongVisAngle = Math.min(90, prevDegDist * (height / currPixDist)) 
            longVisAngle = newLongVisAngle // TODO global
            render(orientQuat)
        } else {
            this.replaceEvent(ev)
        }
    }

    setHandlers(el) {
        el.onpointerdown = (ev) => this.onDownHandler(ev)
        el.onpointermove = (ev) => this.onMoveHandler(ev)
        el.onpointerup = el.onpointercancel = el.onpointerout = el.onpointerleave = (ev) => this.onUpHandler(ev);
    }
}









const isIOS = () => /(iPad|iPhone)/g.test(navigator.userAgent)



const objects = ['Sun', 'Moon', 'Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
const promises = []

promises.push(fetch('./data/stars_vis.json').then(response => response.json()))
promises.push(fetch('./data/planets.json').then(response => response.json()))
for (const obj of objects) {
    promises.push(loadImage('./images/icons/' + obj + '.png'))
}


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



function render(orientQuat) {
    const inverseOrientQuat = Quaternions.inverse(orientQuat)
    const jd = toJd(Date.now())
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    for (const star of stars) {
        let { name, mag, ra, dec } = star
        if (mag < minVisibleMag) {
            // celestial ra from hvg dataset is in hours, so times 15 to get degrees
            const coords = toCanvasCoords(jd, ra*15, dec, inverseOrientQuat)
            if (coords !== null) {
                drawStar(name, mag, coords[0], coords[1])
            }
        }
    }
    {
        const moonLoc = eclipticToEquitorial(Moon.eclipLatLong(jd))
        const coords = toCanvasCoords(jd, moonLoc.ra, moonLoc.dec, inverseOrientQuat)
        if (coords !== null) {
            drawMoon(...coords)
        }
    }
    {
        const sunLoc = eclipticToEquitorial(sunEclipLatLong(jd, earth))
        const coords = toCanvasCoords(jd, sunLoc.ra, sunLoc.dec, inverseOrientQuat)
        if (coords !== null) {
            drawSun(...coords)
        }
    }
    
    for (const p of planets) {
        if (p.name !== "Earth") {
            const planetLoc = eclipticToEquitorial(planetEclipLatLong(jd, p, earth))
            const coords = toCanvasCoords(jd, planetLoc.ra, planetLoc.dec, inverseOrientQuat)
            if (coords !== null) {
                drawPlanet(p, ...coords)
            }
        }
    }
}

function iOSGetOrientationPerms() {
    document.getElementById("request-perms").style.display = 'none';

    if (typeof DeviceOrientationEvent.requestPermission === 'function') {
      DeviceOrientationEvent.requestPermission(true)
        .then(permissionState => {
          if (permissionState === 'granted') {
            window.addEventListener('deviceorientation', () => {
                console.log(event.absolute, event.alpha, event.beta, event.gamma, event.webkitCompassHeading)
                const q = Quaternions.fromAngles(-event.webkitCompassHeading, event.beta, event.gamma)
                console.log('quat', q[0], q[1], q[2], q[3])
                render(q)
            }, true);
          }
        })
        .catch(console.error);
    } 
}

if (isIOS()) {
    const allowButton = document.getElementById("request-perms")
    allowButton.style.display = 'block';
    allowButton.onclick = iOSGetOrientationPerms;
} else {
    const options = { frequency: 30, referenceFrame: "device" };
    const sensor = new AbsoluteOrientationSensor(options);
    sensor.start();
    sensor.addEventListener("reading", () => {
        orientQuat = Quaternions.toInternalQuat(sensor.quaternion)

        //const decFix = Quaternions.fromAngleAxis(2.8, [0, 0, -1])
        //const q = Quaternions.multiply(orientQuat, decFix)

        render(orientQuat)
     });
    sensor.addEventListener("error", (error) => console.log(error));
}

const pinchZoom = new PinchZoom();
pinchZoom.setHandlers(document.getElementById("canvas"));


