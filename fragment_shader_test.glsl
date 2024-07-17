precision highp float;

#define PI 3.1415926538

varying vec3 intersectionPoint;

uniform int maxTraceLevel;

// the mirrors

const int mirrorsN2 = 4;


// uniform bool cylindricalMirrors;
uniform int mirrorType;
uniform float reflectionCoefficient;

uniform vec3 sphereCentre;
uniform float sphereRadius;
uniform bool showSphere;

// background
uniform sampler2D backgroundTexture;

// the camera's wide aperture
uniform float focusDistance;
uniform int noOfRays;
uniform vec3 apertureXHat;
uniform vec3 apertureYHat;
uniform vec3 viewDirection;
uniform float apertureRadius;
uniform float randomNumbersX[100];
uniform float randomNumbersY[100];
// uniform float apertureRadius;

// uniform float xPlaneX[2];

//vec3 xHat = vec3(1.0, 0.0, 0.0);
//vec3 yHat = vec3(0., 1., 0.);
vec3 zHat = vec3(0., 0., 1.);

// // rotation matrix that rotates 2D vectors by the angle alpha (in radians)
// // from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
// mat2 getRotationMatrix(float alpha) {
// 	float s = sin(alpha);
// 	float c = cos(alpha);
// 	return mat2(c, s, -s, c);
// }

// // rotate the 2D vector v by the angle alpha (in radians)
// // from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
// vec2 rotate(vec2 v, float alpha) {
// 	return getRotationMatrix(alpha) * v;
// }

// propagate the ray starting at position s and with direction d to the plane r.n = n0, 
// provided that plane is in the ray's "forward" direction;
// nHat is the *normalised* normal to the plane,
// n0 is the distance of the plane, in the direction n, from the origin,
// p becomes the point where the ray intersects the plane;
// isForward is set to true or false depending on whether the intersection point is forwards or backwards along the ray
void propagateForwardToPlane(
	inout vec3 s, 
	vec3 d, 
	vec3 nHat,
	float n0,
	inout bool isForward
) {
	// calculate the distance in the n direction from the ray start position to the array
	float deltaN = n0 - dot(s, nHat);

	// is the intersection with the plane in the ray's "forward" direction?
	float dN = dot(d, nHat);
	isForward = (dN*deltaN > 0.0);

	// if the intersection is in the forward direction, advance the ray to the plane
	if(isForward) s += d/dN*deltaN;	// set s to the intersection point with the plane
}

// reflect the ray with start point s and direction d off the plane r.nHat = n0,
// provided the intersection with the plane is in the ray's "forward" direction
void planarMirrorReflect(
	inout vec3 s, 
	vec3 d, 
	vec3 nHat,
	float n0,
	inout bool isForward
) {
	propagateForwardToPlane(s, d, nHat, n0, isForward);

	if(isForward) {
		d -= 2.0*d/length(d)*dot(d, nHat);
	}
}

// Calculate the light-ray direction after transmission through a (spherical or cylindrical)
// ideal thin lens or ideal thin (imaging) mirror, or a phase hologram thereof.
// d is the incident light-ray direction;
// pi is a 2D vector containing the vector PI = I-P, i.e. the vector from the principal point P to the intersection point I;
// TO SIMULATE A CYLINDRICAL LENS with (normalised) optical-power direction opdHat, set pi not to (I-P), but to
// the component of (I-P) in the direction of the optical-power direction, opdHat*dot(I-P, opdHat);
// nHat is the normalised normal to the lens/mirror (a factor -1 doesn't matter);
// op is the optical power;
// reflectionFactor should be +1 for a lens and -1 for a mirror;
// returns the outgoing light-ray direction

void lensOrMirrorDeflect(inout vec3 d, vec3 pi, vec3 nHat, float op, float reflectionFactor, bool phaseHologram) {
	vec3 d1, d1T;
	float d1N;
	if(phaseHologram) {
		// phase hologram
		// normalise d such that its length is 1
		d1 = d/length(d);

		d1N = dot(d1, nHat);	// the n component of d1
		vec3 d1T = d1 - nHat*d1N;	// the transverse (perpendicular to nHat) part of d1

		// transverse components of the outgoing light-ray direction
		vec3 dTOut = d1T - pi*op;

		// from the transverse direction, construct a 3D vector by setting the n component such that the length
		// of the vector is 1
		d = dTOut + nHat*reflectionFactor*sign(d1N)*sqrt(1.0 - dot(dTOut, dTOut));
	} else {
		// ideal thin lens/mirror
		// "normalise" the direction such that the magnitude of the n component is 1
		d1 = d/abs(dot(d, nHat));

		d1N = dot(d1, nHat);	// the n component of d1
		vec3 d1T = d1 - nHat*d1N;	// the transverse (perpendicular to nHat) part of d1

		// the 3D deflected direction comprises the transverse components and a n component of magnitude 1
		// and the same sign as d1N = dot(d, nHat)
		d = d1T - pi*op + nHat*reflectionFactor*sign(d1N);
	} 
}

// Calculate the light-ray direction after transmission through a (spherical or cylindrical) lens or lens hologram.
// d is the incident light-ray direction;
// pixy is a 2D vector containing the transverse (x, y) components of the vector I-P,
// i.e. the vector from the principal point P to the intersection point I;
// TO SIMULATE A CYLINDRICAL LENS, replace (I-P).xy with (I-P).pdHat = pdHat * (I-P).pdHat,
// i.e. the parts of (I-P).xy parallel to the optical-power direction
// f is the focal length;
// returns the outgoing light-ray direction

void lensDeflect(inout vec3 d, vec2 pixy, float f, bool idealLens) {
	if(idealLens) {
		// ideal thin lens
		// "normalise" the direction such that the magnitude of the z component is 1
		vec3 d1 = d/abs(d.z);

		// the 3D deflected direction comprises the transverse components and a z component of magnitude 1
		// and the same sign as d.z
		d = vec3(d1.xy - pixy/f, d1.z);
	} else {
		// lens hologram
		// normalise d
		vec3 dN = d/length(d);
		// transverse components of the outgoing light-ray direction
		vec2 dxy = dN.xy - pixy/f;

		// from the transverse direction, construct a 3D vector by setting the z component such that the length
		// of the vector is 1
		d = vec3(dxy, sign(d.z)*sqrt(1.0 - dot(dxy, dxy)));
	}
}

// Pass the current ray (start point s, direction d, brightness factor b) through (or around) a lens.
// The (ideal thin) lens, of focal length f, is in a z plane through centreOfLens.
// It is circular, with the given radius, centred on centreOfLenss.
void passThroughZLens(
	inout vec3 s, 
	inout vec3 d, 
	inout vec4 b,
	vec3 centreOfLens, 
	float radius,
	float focalPower,
	bool idealLens
) {
	bool isForward;
	propagateForwardToPlane(s, d, zHat, centreOfLens.z, isForward);

	if(isForward) {
		// there is an intersection with the plane of this lens in the ray's forward direction

		// does the intersection point lie within the radius?
		vec3 pi = s - centreOfLens;	// vector from the princpal point to the intersection point, p
		float r2 = dot(pi, pi);
		// vec2 pixy = s.xy - centreOfLens.xy;
		// float r2 = dot(pixy, pixy);
		if(r2 < radius*radius) {
			// the intersection point lies inside the radius, so the lens does something to the ray

			// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
			// lensDeflect(d, pixy, focalLength, idealLens);
			lensOrMirrorDeflect(d, pi, 
				vec3(0., 0., 1.),	// nHat 
				// 0.,	// n0
				focalPower,	// f
				-1.,	// reflectionFactor, +1 = transmissive
				false	// phaseHologram
			);

			// lower the brightness factor, giving the light a blue tinge
			b *= vec4(0.9, 0.9, 0.99, 1);
		} 
	}
}

vec4 getColorOfBackground(
	vec3 d
) {
	float l = length(d);
	float phi = atan(d.z, d.x) + PI;
	float theta = acos(d.y/l);
	return texture2D(backgroundTexture, vec2(mod(phi/(2.*PI), 1.0), 1.-theta/PI));
}

bool findNearestIntersectionWithCylinder(	
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// cylinder centre
	float r,	// cylinder radius
	float y_min, // cylinder height
	float y_max, // cylinder height 
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal) {
		
	// for maths see geometry.pdf
	vec3 v = s - c;
	float A = dot(d, d) - d.y*d.y; // we only care about x,z components, subtracting the y part 
	float B = 2.*dot(d, v) - 2.*d.y*v.y;
	float C = dot(v, v) - r*r - v.y*v.y;

	// determinant
	float D = B*B - 4.*A*C;
	
	// is the determinant > 0?
	if (D<0.) {
		// no -- no solution
		return false;
	}

	float sd = sqrt(D);

	// the two solutions of the quadratic equation for delta
	float delta1 = (-B - sd)/(2.*A);
	float delta2 = (-B + sd)/(2.*A);
	
	// the y coordinates of the intersection points corresponding to these solutions
	float y1 = s.y + d.y * delta1; 
	float y2 = s.y + d.y * delta2;

	// calculate the delta that corresponds to the closer intersection in the forward direction
	float DELTA = 1e20;
	if (delta1 > 0. && y1 > y_min && y1 < y_max && delta1 < DELTA) DELTA=delta1;
	
	if (delta2 > 0. && y2 > y_min && y2 < y_max && delta2 < DELTA) DELTA=delta2;

	if (DELTA == 1e20) return false;
	

	intersectionPosition = s + DELTA*d;
	intersectionNormal = intersectionPosition - c;
	intersectionNormal.y = 0.;
	intersectionNormal = normalize(intersectionNormal);
	intersectionDistance = DELTA*length(d);
	return true;
 }

// void phaseHol(vec d, ) {}

bool findNearestIntersectionWithSphere(
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// sphere centre
	float r,	// sphere radius
	out vec3 intersectionPosition,
	out float intersectionDistance
) {
	// for maths see geometry.pdf
	vec3 v = s - c;
	float A = dot(d, d);
	float B = 2.*dot(d, v);
	float C = dot(v, v);

	// calculate the discriminant
	float D= B*B - 4.*A*C;

	if(D < 0.) {
		// the discriminant is negative -- all solutions are imaginary, so there is no intersection
		return false;
	}

	// there is at least one intersection, but is at least one in the forward direction?

	// calculate the square root of the discriminant
	float sd = sqrt(D);

	// try the "-" solution first, as this will be closer, provided it is positive (i.e. in the forward direction)
	float delta = (-B - sd)/(2.*A);
	bool intersection;
	if(delta < 0.) {
		// the delta for the "-" solution is negative, so that is a "backwards" solution; try the "+" solution
		delta = (-B + sd)/(2.*A);

		if(delta < 0.)
			// the "+" solution is also in the backwards direction
			return false;
	}

	// there is an intersection in the forward direction, at
	intersectionPosition = s + delta*d;
	intersectionDistance = delta*length(d);
	return true;
}

// find the closest intersection in the ray's forward direction with either the x, y or z planes
// or any other objects (such as a red sphere)
// s: ray start point (will not be altered)
// d: ray direction
// intersectionPosition: initial value ignored; becomes the position of the intersection
// intersectionDistance: initial value ignored; becomes the distance to the closest intersection point
// objectSetIndex: 0/1/2 if the intersection is with the x/y/z planes, 3 if it is with coloured spheres
// objectIndex: initial value ignored; becomes the index of the object within the object set being intersected
// returns true if an intersection has been found
bool findNearestIntersectionWithObjects(
	vec3 s, 
	vec3 d,
	in int startIntersectionObjectIndex,
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal,
	out int objectIndex
) {
	intersectionDistance = 1e20;	// this means there is no intersection, so far

	// create space for the current...
	vec3 ip;	// ... intersection point, ...
	float id;	// ... intersection distance, ...
	vec3 iN; //... intersection normal...

	objectIndex = startIntersectionObjectIndex;

	// is there an intersection with the cylinder1
	if( showSphere && findNearestIntersectionWithCylinder(s, d, sphereCentre, sphereRadius, -0.2, 0.2, ip, id, iN) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
	if(id < intersectionDistance && (id > 1e-2 || objectIndex != 0)) {
			// the intersection with the z mirrors is the closest one so far
			intersectionPosition = ip;
			intersectionDistance = id;
			intersectionNormal = iN;
			objectIndex = 0;	// sphere
		}
	}
	// is there an intersection with the cylinder2
	if( showSphere && findNearestIntersectionWithCylinder(s, d, sphereCentre, sphereRadius*2., -0.2, 0.2, ip, id, iN) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
		if(id < intersectionDistance && (id > 1e-2 || objectIndex != 1)) {
			// the intersection with the z mirrors is the closest one so far
			intersectionPosition = ip;
			intersectionDistance = id;
			intersectionNormal = iN;
			objectIndex = 1;	// sphere
		}
	}
	
	return (intersectionDistance < 1e20);
}

void main() {
	// first calculate the focusPosition, i.e. the point this pixel is focussed on
	vec3 pv = intersectionPoint - cameraPosition;	// the "pixel view direction", i.e. a vector from the centre of the camera aperture to the point on the object the shader is currently "shading"
	vec3 focusPosition = cameraPosition + focusDistance/abs(dot(pv, viewDirection))*pv;	// see Johannes's lab book 30/4/24 p.174

	// trace <noOfRays> rays
	gl_FragColor = vec4(0, 0, 0, 0);
	vec4 color;
	for(int i=0; i<noOfRays; i++) {
		// the current ray start position, a random point on the camera's circular aperture
		vec3 s = cameraPosition + apertureRadius*randomNumbersX[i]*apertureXHat + apertureRadius*randomNumbersY[i]*apertureYHat;

		// first calculate the current light-ray direction:
		// the ray first passes through focusPosition and then p,
		// so the "backwards" ray direction from the camera to the intersection point is
		vec3 d = focusPosition - s;
		// we normalise this here such that ???
		// d = pv.z/d.z*d;

		// current brightness factor; this will multiply the colour at the end
		vec4 b = vec4(1.0, 1.0, 1.0, 1.0);

		vec3 ip;
		float id;
		vec3 iN;
		float mop;
		vec3 mp;
		vec3 mNHat;
		int oi = -1;
		// int si = -1;
		int tl = maxTraceLevel;	// max trace level
		while(
			(tl-- > 0) &&
			findNearestIntersectionWithObjects(s, d, 
				-1,
				ip,	// out vec3 intersectionPosition
				id,	// out float intersectionDistance
				iN,
				oi	// out int objectIndex
			)
		) {
			s=ip;
			if(oi == 0) { 
				// the first cylinder
				
				b *= vec4(0.6, 1., 1., 1.);
				
				//tl = -10;
			}
			else if(oi == 1) { 
				// the second cylinder
				b *= vec4(0.5882, 0.3176, 0.6627, 1.0);
				//tl = -10;
			}
			// si=oi;
		}
		
		
		if(tl > 0) {
			//color = getColorOfBackground(d);
			color = vec4(1.0);
		} else if(tl != -11) {
			// max number of bounces, but not the sphere
			color = vec4(0.87, 1.0, 0.0, 0.93);
		}

		// finally, multiply by the brightness factor and add to gl_FragColor
		gl_FragColor += b*color;
	}
		
	gl_FragColor /= float(noOfRays);
}