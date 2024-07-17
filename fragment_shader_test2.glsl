precision highp float;

#define PI 3.1415926538

varying vec3 intersectionPoint;

uniform int maxTraceLevel;

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

bool findNearestIntersectionWithCylinder(	
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// cylinder centre
	float r,	// cylinder radius
	float y_min, // cylinder height
	float y_max, // cylinder height 
	bool startPointIsIntersectionWithCylinder,
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal
) {
		
	// for maths see geometry.pdf
	vec2 v2 = s.xz - c.xz;
	vec2 d2 = d.xz;
	float A = dot(d2, d2);
	float B = 2.*dot(d2, v2);

	float delta = 1e20;

	if(startPointIsIntersectionWithCylinder) {
		// if the ray start point lies son the cylinder, then C = 0
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

	intersectingObjectIndex = startIntersectionObjectIndex;

	// is there an intersection with the cylinder1
	if( true && findNearestIntersectionWithCylinder(s, d, sphereCentre, sphereRadius, -0.2, 0.2, intersectingObjectIndex == 0, intersectionPosition, intersectionDistance, intersectionNormal) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
		if(intersectionDistance < closestIntersectionDistance && (intersectionDistance > 1e-5 || intersectingObjectIndex != 0)) {
			// the intersection with the z mirrors is the closest one so far
			closestIntersectionPosition = intersectionPosition;
			closestIntersectionDistance = intersectionDistance;
			closestIntersectionNormal = intersectionNormal;
			intersectingObjectIndex = 0;	// sphere
		}
	}
	// is there an intersection with the cylinder2
	if( true && findNearestIntersectionWithCylinder(s, d, sphereCentre, 2.*sphereRadius, -0.1, .1, intersectingObjectIndex == 1, intersectionPosition, intersectionDistance, intersectionNormal) ) {
		// yes, there is an intersection with the cylinder
		// if there either no intersection already, or, if there is one, is it closer than the closest intersection so far?
		if(intersectionDistance < closestIntersectionDistance && (intersectionDistance > 1e-5 || intersectingObjectIndex != 1)) {
			// the intersection with the z mirrors is the closest one so far
			closestIntersectionPosition = intersectionPosition;
			closestIntersectionDistance = intersectionDistance;
			closestIntersectionNormal = intersectionNormal;
			intersectingObjectIndex = 1;	// sphere
		}
	}
	
	return (closestIntersectionDistance < 1e20);
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
				oi,
				ip,	// out vec3 intersectionPosition
				id,	// out float intersectionDistance
				iN,
				oi	// out int objectIndex
			)
		) {
			if(oi == 0) { 
				// the first cylinder
				b *= vec4(1., .4, .4, 1.);
			}
			else if(oi == 1) { 
				// the second cylinder
				b *= vec4(.4, 0.4, 1., 1.0);
			}

			s=ip;
		}
		
		// if(tl > 0) {
		// 	color = getColorOfBackground(d);
		// } 
		// else if(tl != -11) {
		// 	// max number of bounces
		// 	color = vec4(0.0, 0.0, 0.0, 1.0);
		// }
		color = vec4(1., 1., 1., 1.);

		// finally, multiply by the brightness factor and add to gl_FragColor
		gl_FragColor += b*color;
	}
		
	gl_FragColor /= float(noOfRays);
}