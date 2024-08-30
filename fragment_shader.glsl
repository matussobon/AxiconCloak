precision highp float;

#define PI 3.1415926538

#define NONE -1

// object types
#define OBJECT_TYPE_RECTANGLE 0

// surface types
#define SURFACE_TYPE_LENS 0
#define SURFACE_TYPE_COLOR 1

// lens types
#define LENS_TYPE_IDEAL 0
#define LENS_TYPE_HOLOGRAM 1


varying vec3 intersectionPoint;

uniform int maxTraceLevel;

uniform float rotAngle;
uniform bool showLens;
uniform float cornerDistance; // distance between the lenses

uniform float sphereRadius;
uniform bool showSphere;
uniform float sphereHeight;
uniform vec3 sphereCentre;
// show/hide the whole Axicon Cloak
uniform bool showCloak;


// Axicon Cloak centre 
uniform bool cloakCentre;
uniform float yShift;
uniform float phaseShift;

// to do: show/hide the individual cylinders
uniform bool showOuterCylinder;
uniform bool showInnerCylinder;

// outer cylinder properties
uniform float outerRadius; 
uniform float outerHeightNegative;
uniform float outerHeightPositive;
uniform float outerYcoord;


// inner cylinder properties
uniform float innerRadius;
uniform float innerHeightNegative; 
uniform float innerHeightPositive;
uniform float innerYcoord;


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

// camera aperture struct
struct CameraAperture {
	vec3 viewDirection;
	vec3 apertureXHat;
	vec3 apertureYHat;
	float focusDistance;
	float apertureRadius;
	float randomNumbersX[100];
	float randomNumbersY[100];
	int noOfRays;
};

uniform CameraAperture Camera;

struct Ray {
	vec3 origin;
	vec3 direction;
};

struct addObject {
	bool visible;
	vec3 centre;
	float size; 
};
// size variable has different meaning for each object i.e.
// Sphere, Cylinder: size = radius 
// Rectangle: size = side-length

uniform addObject Sphere;
uniform addObject Rectangle;
uniform addObject Cylinder;

vec3 uSpan = vec3(1,0,0);
vec3 vSpan = vec3(0,1,0);

struct Rectangles { 
	bool visible;
	vec3 corner;
    vec3 uSpanVector;
    vec3 vSpanVector;
	float uSize;
	float vSize;
    int surfaceType;
    int surfaceIndex; 
};
uniform Rectangles rectangles[ 3 ];


struct LensSurface {
	vec3 principalPoint;
    vec3 opticalAxisDirection;
    float focalLength;
    float transmissionCoefficient;
	int lensType;
};
uniform LensSurface lensSurfaces[ 3 ];


struct Colour {
	vec4 color;
};
uniform Colour colors[ 2 ];

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

bool findNearestIntersectionWithLens(
	vec3 s, // ray start point, origin 
	vec3 d, // ray direction 
	vec3 lensCorner,
	float lensSize, // assuming rectangular aperture
	float rotAngle, // for now rotation around Y-axis 
	out vec3 intersectionPosition,
	out float intersectionDistance
) {
	float cosTheta = cos(radians(-rotAngle));
	float sinTheta = sin(radians(-rotAngle));

	mat3 rotMatrix = mat3(
		cosTheta, 0., sinTheta,
		0., 	 1., 	0.,
		-sinTheta,0., cosTheta
	);

	vec3 lensNormal = rotMatrix * vec3 (0, 0, 1);

	// calculate the span vectors
	vec3 uSpan = rotMatrix * vec3(1., 0., 0.);
	vec3 vSpan = rotMatrix * vec3(0., 1., 0.);

	// if the ray is parallel to the lens surface there is no intersection 
	if (dot(d, lensNormal)==0.) {
		return false;
	}
	// calculate delta to check for intersections 
	float delta = dot(lensCorner - s, lensNormal)/(dot(d, lensNormal));
	intersectionPosition = s + delta*d;

	if (delta<0.) {
		return false;
	} 
	float uProj = dot((intersectionPosition - lensCorner),uSpan);
	float vProj = dot(intersectionPosition - lensCorner, vSpan);

	if (uProj<0. || uProj>lensSize || vProj<0. || vProj>lensSize){
		return false;
	}

	intersectionDistance = delta*length(d);
	return true;
}

bool findNearestIntersectionWithRectangle(
	vec3 s, // ray start point, origin 
	vec3 d, // ray direction 
	vec3 corner,
	vec3 uSpanVector,	// normalised!
	vec3 vSpanVector,	// normalised!
	float uSize,
	float vSize,
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal
) {
	vec3 normal = normalize(cross(uSpanVector, vSpanVector));

	// if the ray is parallel to the lens surface there is no intersection 
	if (dot(d, normal)==0.) {
		return false;
	}

	// calculate delta to check for intersections 
	float delta = dot(corner - s, normal)/(dot(d, normal));
	intersectionPosition = s + delta*d;

	if (delta<0.) {
		return false;
	} 

	vec3 ci = intersectionPosition - corner;
	float uProj = dot(ci, uSpanVector);
	float vProj = dot(ci, vSpanVector);

	if (uProj<0. || uProj>uSize || vProj<0. || vProj>vSize){
		return false;
	}

	intersectionDistance = delta*length(d);
	intersectionNormal = normal;
	return true;
}

bool findNearestIntersectionWithCylinder(	
	vec3 s, 	// ray start point
	vec3 d, 	// ray direction
	vec3 c,		// cylinder centre
	float r,	// cylinder radius
	float yShift,
	float yCoord,
	float y_min, // cylinder height
	float y_max, // cylinder height 
	bool startPointIsIntersectionWithCylinder,
	out vec3 intersectionPosition,
	out float intersectionDistance,
	out vec3 intersectionNormal
) {
	y_min = y_min + yShift + yCoord;
	y_max = y_max + yShift+ yCoord;
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
	in int startObjectType,	// rectangle = 0
	in int startObjectIndex,
	out vec3 closestIntersectionPosition,
	out float closestIntersectionDistance,
	out vec3 closestIntersectionNormal,
	out int closestIntersectionObjectType,
	out int closestIntersectionObjectIndex,
	out int closestIntersectionSurfaceType,
	out int closestIntersectionSurfaceIndex
) {
	closestIntersectionDistance = 1e20;	// this means there is no intersection, so far

	// create space for the current...
	vec3 intersectionPosition;	// ... intersection point, ...
	float intersectionDistance;	// ... intersection distance, ...
	vec3 intersectionNormal; //... intersection normal...

	// is there an intersection with the array of rectangles?
	int i = 0;
	// go through all the rectangles
	while( i < rectangles.length() ) {
		// check if the ray started on this rectangle
		if( ((startObjectType != OBJECT_TYPE_RECTANGLE) || (startObjectIndex != i)) && rectangles[i].visible ) {
			// check if the intersection with it is closer than the closest one so far
			if (findNearestIntersectionWithRectangle(
				s, d, rectangles[i].corner, rectangles[i].uSpanVector, rectangles[i].vSpanVector, 
				rectangles[i].uSize, rectangles[i].vSize, intersectionPosition, intersectionDistance, intersectionNormal)
			) {
				if(intersectionDistance < closestIntersectionDistance) {
					closestIntersectionPosition = intersectionPosition;
					closestIntersectionDistance = intersectionDistance;
					closestIntersectionNormal = intersectionNormal;
					closestIntersectionObjectType = OBJECT_TYPE_RECTANGLE;	// rectangle
					closestIntersectionObjectIndex = i;
					closestIntersectionSurfaceType = rectangles[i].surfaceType;
					closestIntersectionSurfaceIndex = rectangles[i].surfaceIndex;
				}
			}
		}
		i++;
	}

	// now do the same for all other object types

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

// return the light-ray direction after transmission through a lens of focal length f with 
// (normalised) optical-axis direction axisHat;
// d is the initial light-ray direction;
// p2i is a vector from the lens's principal point to the intersection point (and this should be perpendicular to axisHat)
vec3 lensDeflect(vec3 d, vec3 p2i, vec3 axisHat, float f, bool idealLens) {

    if(idealLens) {
        // ideal thin lens
		// d' \propto d - d_axial / f (I - P) -- see geometry.pdf
		return d - (dot(d, axisHat)/f)*p2i;
	}
}

void main() {
	Ray LightRay;
	// first calculate the focusPosition, i.e. the point this pixel is focussed on
	vec3 pv = intersectionPoint - cameraPosition;	// the "pixel view direction", i.e. a vector from the centre of the camera aperture to the point on the object the shader is currently "shading"
	vec3 focusPosition = cameraPosition + focusDistance/abs(dot(pv, viewDirection))*pv;	// see Johannes's lab book 30/4/24 p.174
	
	// trace <noOfRays> rays
	gl_FragColor = vec4(0, 0, 0, 0);
	vec4 color;
	for(int i=0; i<noOfRays; i++) {
		// the current ray start position, a random point on the camera's circular aperture

		LightRay.origin = cameraPosition + Camera.apertureRadius*Camera.randomNumbersX[i]*Camera.apertureXHat + Camera.apertureRadius*Camera.randomNumbersY[i]*Camera.apertureYHat;

		// LightRay.origin = cameraPosition + apertureRadius*randomNumbersX[i]*apertureXHat + apertureRadius*randomNumbersY[i]*apertureYHat;

		// first calculate the current light-ray direction:
		// the ray first passes through focusPosition and then p,
		// so the "backwards" ray direction from the camera to the intersection point is
		LightRay.direction = focusPosition - LightRay.origin;
		// we normalise this here such that ???
		// d = pv.z/d.z*d;

		// current brightness factor; this will multiply the colour at the end
		vec4 b = vec4(1.0, 1.0, 1.0, 1.0);

		// int si = -1;
		int traceLevel = maxTraceLevel;	// max trace level

		vec3 intersectionPosition;
		float intersectionDistance;
		vec3 intersectionNormal;
		int intersectionObjectType = NONE;
		int intersectionObjectIndex = NONE;
		int intersectionSurfaceType;
		int intersectionSurfaceIndex;

		while(
			(traceLevel-- > 0) &&
			findNearestIntersectionWithObjects(LightRay.origin, LightRay.direction, 
				 intersectionObjectType,
				 intersectionObjectIndex,
				 intersectionPosition,
				 intersectionDistance,
				 intersectionNormal,
				 intersectionObjectType,
				 intersectionObjectIndex,
				 intersectionSurfaceType,
				 intersectionSurfaceIndex
			)
		) {
			// there is an intersection; what type of surface is it?
			if(intersectionSurfaceType == SURFACE_TYPE_LENS) {
				// if(lensSurfaces[intersectionSurfaceIndex].lensType == LENS_TYPE_IDEAL) {
					LightRay.direction = lensDeflect(
						LightRay.direction, 
						intersectionPosition-lensSurfaces[intersectionSurfaceIndex].principalPoint, 
						lensSurfaces[intersectionSurfaceIndex].opticalAxisDirection, 
						lensSurfaces[intersectionSurfaceIndex].focalLength, 
						lensSurfaces[intersectionSurfaceIndex].lensType == LENS_TYPE_IDEAL
					);
					b.rgb *= lensSurfaces[intersectionSurfaceIndex].transmissionCoefficient;
				// } else {
				// 	LightRay.direction = phaseHologram(
				// 		LightRay.direction, 
				// 		lensSurfaces[intersectionSurfaceIndex].opticalAxisDirection,
				// 		0.5
				// 	);
				// }
				LightRay.origin=intersectionPosition;
			} 
			else if(intersectionSurfaceType == SURFACE_TYPE_COLOR) {
				color = colors[intersectionSurfaceIndex].color;
				traceLevel = -10;
			}
		}
		
		if(traceLevel > 0) {
			color = getColorOfBackground(LightRay.direction);
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