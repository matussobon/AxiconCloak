precision highp float;

#define PI 3.1415926538

varying vec3 intersectionPoint;

uniform int maxTraceLevel;


uniform float sphereRadius;
uniform bool showSphere;
uniform float sphereHeight;
uniform vec3 sphereCentre;
// show/hide the whole Axicon Cloak
uniform bool showCloak;

// Axicon Cloak centre 
uniform bool cloakCentre;
uniform float yShift;

// to do: show/hide the individual cylinders
uniform bool showOuterCylinder;
uniform bool showInnerCylinder;

// outer cylinder properties
uniform float outerRadius; 
uniform float outerHeightNegative;
uniform float outerHeightPositive;
uniform float outerPhaseShift;

// inner cylinder properties
uniform float innerRadius;
uniform float innerHeightNegative; 
uniform float innerHeightPositive;
uniform float innerPhaseShift;

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


vec3 zHat = vec3(0., 0., 1.);

vec4 getColorOfBackground(
	vec3 d
) {
	float l = length(d);
	float phi = atan(d.z, d.x) + PI;
	float theta = acos(d.y/l);
	return texture2D(backgroundTexture, vec2(mod(phi/(2.*PI), 1.0), 1.-theta/PI));
}

bool findNearestIntersectionWithSphere(
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// sphere centre
	float y,  	// y coordinate of sphere centre 
	float yShift,
	float r,	// sphere radius
	out vec3 intersectionPosition,
	out float intersectionDistance
) {
	c = c + vec3(0, y+yShift , 0);
	// for maths see geometry.pdf
	vec3 v = s - c;
	float A = dot(d, d);
	float B = 2.*dot(d, v);
	float C = dot(v, v) - r*r;

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
	//bool intersection;
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

bool findNearestIntersectionWithCylinder(	
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// cylinder centre
	float r,	// cylinder radius
	float yShift,
	float y_min, // cylinder height
	float y_max, // cylinder height 
	bool startPointIsIntersectionWithCylinder,
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal
) {
	y_min = y_min + yShift;
	y_max = y_max + yShift;
	// for maths see geometry.pdf
	vec2 v2 = s.xz - c.xz;
	vec2 d2 = d.xz;
	float A = dot(d2, d2);
	float B = 2.*dot(d2, v2);

	float delta = 1e20;

	if(startPointIsIntersectionWithCylinder) {
		// if the ray start point lies on the cylinder, then C = 0
		delta = -B/A;
		if(delta > 0.) {
			float y = s.y + delta*d.y;
			if(y_min <= y && y <= y_max) {
				intersectionPosition = s + delta*d;
				intersectionNormal = intersectionPosition - c;
				intersectionNormal.y = 0.;
				intersectionNormal = normalize(intersectionNormal);
				intersectionDistance = delta*length(d);
				return true;
			}
		}
		return false;
	} 
	
	float C = dot(v2, v2) - r*r;

	// determinant
	float D = B*B - 4.*A*C;
	
	// is the determinant > 0?
	if (D<0.) return false; // no -- no solution

	float sd = sqrt(D);

	// the two solutions of the quadratic equation for delta
	float delta1 = (-B - sd)/(2.*A);
	float delta2 = (-B + sd)/(2.*A);
	
	// the y coordinates of the intersection points corresponding to these solutions
	float y1 = s.y + d.y * delta1; 
	float y2 = s.y + d.y * delta2;

	// calculate the delta that corresponds to the closer intersection in the forward direction
	if (delta1 > 0. && y1 > y_min && y1 < y_max && delta1 < delta) delta=delta1;
	
	if (delta2 > 0. && y2 > y_min && y2 < y_max && delta2 < delta) delta=delta2;

	if (delta == 1e20) return false;

	intersectionPosition = s + delta*d;
	intersectionNormal = intersectionPosition - c;
	intersectionNormal.y = 0.;
	intersectionNormal = normalize(intersectionNormal);
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
	out vec3 closestIntersectionPosition,
	out float closestIntersectionDistance,
	out vec3 closestIntersectionNormal,
	out int intersectingObjectIndex
) {
	closestIntersectionDistance = 1e20;	// this means there is no intersection, so far

	// create space for the current...
	vec3 intersectionPosition;	// ... intersection point, ...
	float intersectionDistance;	// ... intersection distance, ...
	vec3 intersectionNormal; //... intersection normal...

	// is there an intersection with the cylinder1
	if( showCloak && showInnerCylinder && findNearestIntersectionWithCylinder(s, d, sphereCentre, innerRadius, yShift, innerHeightNegative, innerHeightPositive, startIntersectionObjectIndex == 0, intersectionPosition, intersectionDistance, intersectionNormal) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
		if(intersectionDistance < closestIntersectionDistance && (intersectionDistance > 1e-5 || startIntersectionObjectIndex != 0)) {
			// the intersection with the z mirrors is the closest one so far
			closestIntersectionPosition = intersectionPosition;
			closestIntersectionDistance = intersectionDistance;
			closestIntersectionNormal = intersectionNormal;
			intersectingObjectIndex = 0;	// sphere
		}
	}
	// is there an intersection with the cylinder2
	if( showCloak && showOuterCylinder && findNearestIntersectionWithCylinder(s, d, sphereCentre, outerRadius,yShift, outerHeightNegative + yShift, outerHeightPositive + yShift, startIntersectionObjectIndex == 1, intersectionPosition, intersectionDistance, intersectionNormal) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
		if(intersectionDistance < closestIntersectionDistance && (intersectionDistance > 1e-5 || startIntersectionObjectIndex != 1)) {
			// the intersection with the z mirrors is the closest one so far
			closestIntersectionPosition = intersectionPosition;
			closestIntersectionDistance = intersectionDistance;
			closestIntersectionNormal = intersectionNormal;
			intersectingObjectIndex = 1;	// sphere
		}
	}

	// if (showSphere && findNearestIntersectionWithSphere(s, d, sphereCentre, sphereRadius, intersectionPosition, intersectionDistance))
	if( showSphere && findNearestIntersectionWithSphere(s, d, sphereCentre, sphereHeight, yShift, sphereRadius, intersectionPosition, intersectionDistance) ) {
		if (intersectionDistance < closestIntersectionDistance)  {
			closestIntersectionPosition = intersectionPosition;
			closestIntersectionDistance = intersectionDistance;
			intersectingObjectIndex = 2;
		}
	}

	return (closestIntersectionDistance < 1e20);
}

// d - incident ray direction 
// closestIntersectionNormal - normal to the surface with of norm = 1
// deltaKy - phase shift of the hologram
// ideally, we just use this function for each cylinder and get the outgoing ray direction
// by setting the deltaKy with opposite sign 
vec3 phaseHologram(vec3 d, vec3 closestIntersectionNormal, float deltaKy) {

	//normalize the the ray direction vector d
	vec3 dNorm = normalize(d);
	//calculate the projection of dNorm onto the closestIntersectionNormal
	float dNormProj = dot(dNorm, closestIntersectionNormal);

	
	vec3 dTransverse = dNorm - closestIntersectionNormal*dNormProj; 
	vec3 dPrimeTransverse = dTransverse +  vec3 (0.0, deltaKy ,0.0); 
	
	//check for evanescence
	float dPrimeTransverseNorm = dot(dPrimeTransverse,dPrimeTransverse);
	if (dPrimeTransverseNorm > 1.) {
		//if true
		return reflect(dNorm, closestIntersectionNormal);
		
	}


	float dPrimeNormal = sqrt(1.0-dPrimeTransverseNorm);

	vec3 dPrime = dPrimeTransverse  + sign(dNormProj)*dPrimeNormal*closestIntersectionNormal;
	return dPrime;
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
		int oi = -1;
		// int si = -1;
		int tl = maxTraceLevel;	// max trace level
		while(
			(tl-- > 0) &&
			findNearestIntersectionWithObjects(s, d, 
				oi,
				ip,	// out vec3 intersectionPosition
				id,	// out float intersectionDistance
				iN,
				oi	// out int objectIndex
			)
		) {
			if(oi == 0) { 
				// the first cylinder
				// b *= vec4(1., .5, .5, 1.);
				d = phaseHologram(d, iN, outerPhaseShift );
			}
			else if(oi == 1) { 
				// the second cylinder
				// b *= vec4(.5, 0.5, 1., 1.0);
				d = phaseHologram(d, iN, innerPhaseShift);
			}

			else if (oi == 2) 
			{
				// the sphere 
				color = vec4(1., 0., 0., 1.);
				//b *= vec4(.5, 0.5, 1., 1.0);
				tl = -10;
			}

			s=ip;
		}
		
		if(tl > 0) {
			color = getColorOfBackground(d);
		} 
		// else if(tl != -11) {
		// // 	// max number of bounces
		// color = vec4(0.0, 0.0, 0.0, 1.0);
		// }
		// color = vec4(1., 1., 1., 1.);

		// finally, multiply by the brightness factor and add to gl_FragColor
		gl_FragColor += b*color;
	}
		
	gl_FragColor /= float(noOfRays);
}