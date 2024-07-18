// This code is based on three.js, which comes with the following license:
//
// The MIT License
//
// Copyright © 2010-2024 three.js authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
import * as THREE from "three";

import { GUI } from "three/addons/libs/lil-gui.module.min.js";
import { DragControls } from "three/addons/controls/DragControls.js";
import { VRButton } from "three/addons/webxr/VRButton.js";

import { HTMLMesh } from "three/addons/interactive/HTMLMesh.js";
import { InteractiveGroup } from "three/addons/interactive/InteractiveGroup.js";
import { XRControllerModelFactory } from "three/addons/webxr/XRControllerModelFactory.js";
// import { createMeshesFromInstancedMesh } from 'three/examples/jsm/utils/SceneUtils.js';

import {
  background2String,
  resonatorType2String,
  mirrorType2String,
  showSphere2String,
  guiMeshVisible2String,
  createInfo,
  refreshInfo,
  toggleInfoVisibility,
  showCloak2String,
  showInnerCylinder2String,
  showOuterCylinder2String,
} from "/InfoFunctions.js";
import {
  postStatus,
  createStatus,
  screenChanged,
  onWindowResize,
} from "/statusFunctions.js";
import { loadBackgroundImage, toggleFullscreen } from "./CosmeticsFunctions.js";
import {
  takePhoto,
  showStoredPhoto,
  showLivePhoto,
  deleteStoredPhoto,
} from "./PhotoFunctions.js";
import { addOrbitControls, pointForward } from "./CameraFunctions.js";

// import vertexShaderCode from "./vertex_shader_test.glsl?raw";
//'./shader.glsl?raw'

// import fragmentShaderCode from "./fragment_shader_test.glsl";
// import fragmentShaderCode from "./fragment_shader_test3.glsl?raw";
//https://www.youtube.com/watch?v=RDughHM9qoE
//https://www.npmjs.com/package/vite-plugin-glsl

let scene;
let renderer;
let backgroundTexture;
let orbitControls;
let dragControls;
let raytracingSphere;

let mirrorsN2 = 4; // max number of mirrors in each array

let xMirrorsN;
let xMirrorsX; // {x[0], x[1], ...}; note that we require x[0] <= x[1] <= x[2] ...!
let xMirrorsY1; // {y1[0], y1[1], ...}
let xMirrorsY2; // {y2[0], y2[1], ...}
let xMirrorsZ1; // {z1[0], z1[1], ...}
let xMirrorsZ2; // {z2[0], z2[1], ...}
let xMirrorsP; // {P[0], P[1], ...} Principal points
let xMirrorsOP; // {op[0], op[1], ...} optical powers

let yMirrorsN;
let yMirrorsY; // {y[0], y[1], ...}; note that we require y[0] <= y[1] <= y[2] ...!
let yMirrorsX1; // {x1[0], x1[1], ...}
let yMirrorsX2; // {x2[0], x2[1], ...}
let yMirrorsZ1; // {z1[0], z1[1], ...}
let yMirrorsZ2; // {z2[0], z2[1], ...}
let yMirrorsP; // {P[0], P[1], ...} Principal points
let yMirrorsOP; // {op[0], op[1], ...} optical powers

let zMirrorsN;
let zMirrorsZ; // {z[0], z[1], ...}; note that we require z[0] <= z[1] <= z[2] ...!
let zMirrorsX1; // {x1[0], x1[1], ...}
let zMirrorsX2; // {x2[0], x2[1], ...}
let zMirrorsY1; // {y1[0], y1[1], ...}
let zMirrorsY2; // {y2[0], y2[1], ...}
let zMirrorsP; // {P[0], P[1], ...} Principal points
let zMirrorsOP; // {op[0], op[1], ...} optical powers

let y1 = -0.5;
let y2 = 0.5;

let sphereRadius = 0.1;
let sphereHeight = 0;

let outerRadius = 0.3;
let outerHeightNegative = -0.1;
let outerHeightPositive = 0.1;

let innerRadius = 0.1;
let innerHeightNegative = -0.2;
let innerHeightPositive = 0.2;

let outerPhaseShift = 0.6;
let innerPhaseShift = -0.6;

let raytracingSphereRadius = 100.0;

let autofocus = false;

// the menu
let gui;
let GUIParams;
let autofocusControl,
  focusDistanceControl,
  resonatorTypeControl,
  mirrorTypeControl,
  opz0Control,
  opz1Control,
  z0Control,
  z1Control,
  resonatorYControl,
  cylindricalMirrorsControl,
  backgroundControl,
  vrControlsVisibleControl,
  showSphereControl;

let showCloakControl;
let showInnerCylinderControl;
let showOuterCylinderControl;

let GUIMesh;
let showGUIMesh;
// let meshRotationX = -Math.PI/4, meshRotationY = 0, meshRotationZ = 0;

// true if stored photo is showing
let storedPhoto;

const infoObject = {
  showingStoredPhoto: false,
  storedPhotoInfoString: undefined,
  resonatorType: 1, // 0 = single canonical resonator in x direction, 1 = crossed canonical resonators in x and z directions, 2 = Penrose cavity
  raytracingSphereShaderMaterial: undefined,
  background: 0,
  camera: undefined,
  fovScreen: 68,
  x1: -0.5,
  x2: 0.5,
  z1: -0.5,
  z2: 0.5,
  resonatorY: 0.0, // lift the resonator up to eye level (in case of VR only)
  xMirrorX1OP: 0,
  xMirrorX2OP: 0,
  zMirrorZ1OP: 0,
  zMirrorZ2OP: 0,
  sphereCentre: new THREE.Vector3(0, 0, 0),
  apertureRadius: 0.0, // camera with wide aperture
  atanFocusDistance: Math.atan(3e8),
  noOfRays: 1,
  storedPhotoDescription: undefined,
};

init();
animate();

function init() {
  // create the info element first so that any problems can be communicated
  createStatus();

  scene = new THREE.Scene();
  // scene.background = new THREE.Color( 'skyblue' );
  let windowAspectRatio = window.innerWidth / window.innerHeight;
  infoObject.camera = new THREE.PerspectiveCamera(
    infoObject.fovScreen,
    windowAspectRatio,
    0.1,
    2 * raytracingSphereRadius + 1
  );
  infoObject.camera.position.z = 1;
  screenChanged(renderer, infoObject.camera, infoObject.fovScreen);

  renderer = new THREE.WebGLRenderer({
    antialias: true,
    preserveDrawingBuffer: true,
  });
  renderer.setPixelRatio(window.devicePixelRatio);
  renderer.setSize(window.innerWidth, window.innerHeight);
  renderer.xr.enabled = true;
  document.body.appendChild(VRButton.createButton(renderer)); // for VR content
  document.body.appendChild(renderer.domElement);
  // document.getElementById('livePhoto').appendChild( renderer.domElement );

  backgroundTexture = loadBackgroundImage(
    infoObject.background,
    backgroundTexture
  );

  initMirrors();
  addRaytracingSphere();

  // user interface

  addEventListenersEtc();
  orbitControls = addOrbitControls(infoObject.camera, renderer.domElement); //instead of addOrbitControls();

  // the controls menu
  // refreshGUI();
  createGUI();

  // addDragControls();

  // check if VR is supported (see https://developer.mozilla.org/en-US/docs/Web/API/XRSystem/isSessionSupported)...
  // if (navigator.xr) {
  if ("xr" in navigator) {
    // renderer.xr.enabled = false;
    // navigator.xr.isSessionSupported("immersive-vr").then((isSupported) => {
    navigator.xr.isSessionSupported("immersive-vr").then(function (supported) {
      if (supported) {
        // ... and enable the relevant features
        renderer.xr.enabled = true;
        // use renderer.xr.isPresenting to find out if we are in XR mode -- see https://threejs.org/docs/#api/en/renderers/webxr/WebXRManager
        // (and https://threejs.org/docs/#api/en/renderers/WebGLRenderer.xr, which states that rendereor.xr points to the WebXRManager)
        document.body.appendChild(VRButton.createButton(renderer)); // for VR content
        addXRInteractivity();
      }
    });
  }

  createInfo();

  refreshInfo(infoObject);
}

function animate() {
  renderer.setAnimationLoop(render);
}

function render() {
  // requestAnimationFrame( animate );

  // stats.begin();

  if (!infoObject.showingStoredPhoto) {
    // update uniforms
    updateUniforms();
    renderer.render(scene, infoObject.camera);
  }

  // stats.end();
}

function updateUniforms() {
  // // are we in VR mode?
  let deltaY;
  // if(renderer.xr.enabled && renderer.xr.isPresenting) {
  deltaY = infoObject.resonatorY;
  // } else {
  // 	deltaY = 0;
  // }

  GUIMesh.position.y = deltaY - 1;

  switch (infoObject.resonatorType) {
    case 0: // 0 = no resonator:
      xMirrorsN = 0;
      yMirrorsN = 0;
      zMirrorsN = 0;
      break;
    case 2: // 2 = crossed canonical resonators in x and z directions
      xMirrorsN = 2;
      yMirrorsN = 0;
      zMirrorsN = 2;

      for (let i = 0; i < xMirrorsN; i++) {
        xMirrorsY1[i] = y1 + deltaY;
        xMirrorsY2[i] = y2 + deltaY;
        xMirrorsZ1[i] = infoObject.z1;
        xMirrorsZ2[i] = infoObject.z2;
        xMirrorsP[i].y = deltaY;
      }
      xMirrorsX[0] = infoObject.x1;
      xMirrorsP[0].x = infoObject.x1;
      xMirrorsOP[0] = infoObject.xMirrorX1OP;
      xMirrorsX[1] = infoObject.x2;
      xMirrorsP[1].x = infoObject.x2;
      xMirrorsOP[1] = infoObject.xMirrorX2OP;

      for (let i = 0; i < zMirrorsN; i++) {
        zMirrorsY1[i] = y1 + deltaY;
        zMirrorsY2[i] = y2 + deltaY;
        zMirrorsX1[i] = infoObject.x1;
        zMirrorsX2[i] = infoObject.x2;
        zMirrorsP[i].y = deltaY;
      }
      zMirrorsZ[0] = infoObject.z1;
      zMirrorsP[0].z = infoObject.z1;
      zMirrorsOP[0] = infoObject.zMirrorZ1OP;
      zMirrorsZ[1] = infoObject.z2;
      zMirrorsP[1].z = infoObject.z2;
      zMirrorsOP[1] = infoObject.zMirrorZ2OP;
      break;
    case 3: // 3 = Penrose cavity
      xMirrorsN = 4;
      yMirrorsN = 0;
      zMirrorsN = 4;

      let z3 = infoObject.z1 - 2 / Math.abs(infoObject.zMirrorZ1OP);
      let z4 = infoObject.z2 + 2 / Math.abs(infoObject.zMirrorZ2OP);

      // the x mirrors
      for (let i = 0; i < xMirrorsN; i++) {
        xMirrorsY1[i] = y1 + deltaY;
        xMirrorsY2[i] = y2 + deltaY;
        xMirrorsP[i].y = deltaY;
      }

      xMirrorsZ1[0] = z3;
      xMirrorsZ2[0] = z4;
      xMirrorsX[0] = 2 * infoObject.x1;
      xMirrorsP[0].x = 2 * infoObject.x1;
      xMirrorsOP[0] = 0; // the outer left mirror
      xMirrorsZ1[1] = infoObject.z1;
      xMirrorsZ2[1] = infoObject.z2;
      xMirrorsX[1] = infoObject.x1;
      xMirrorsP[1].x = infoObject.x1;
      xMirrorsOP[1] = 0; // the inner left mirror
      xMirrorsZ1[2] = infoObject.z1;
      xMirrorsZ2[2] = infoObject.z2;
      xMirrorsX[2] = -infoObject.x1;
      xMirrorsP[2].x = -infoObject.x1;
      xMirrorsOP[1] = 0; // the inner right mirror
      xMirrorsZ1[3] = z3;
      xMirrorsZ2[3] = z4;
      xMirrorsX[3] = -2 * infoObject.x1;
      xMirrorsP[3].x = -2 * infoObject.x1;
      xMirrorsOP[3] = 0; // the outer right mirror

      // the z mirrors
      for (let i = 0; i < zMirrorsN; i++) {
        zMirrorsY1[i] = y1 + deltaY;
        zMirrorsY2[i] = y2 + deltaY;
        zMirrorsP[i].y = deltaY;
      }

      zMirrorsX1[0] = 2 * infoObject.x1;
      zMirrorsX2[0] = -2 * infoObject.x1;
      zMirrorsZ[0] = z3;
      zMirrorsP[0].z = z3;
      zMirrorsOP[0] = infoObject.zMirrorZ1OP; // the outer top mirror
      zMirrorsX1[1] = 2 * infoObject.x1;
      zMirrorsX2[1] = infoObject.x1;
      zMirrorsZ[1] = 0;
      zMirrorsP[1].z = 0;
      zMirrorsOP[1] = 0; // the inner top mirror
      zMirrorsX1[2] = -infoObject.x1;
      zMirrorsX2[2] = -2 * infoObject.x1;
      zMirrorsZ[2] = 0;
      zMirrorsP[2].z = 0;
      zMirrorsOP[2] = 0; // the inner bottom mirror
      zMirrorsX1[3] = 2 * infoObject.x1;
      zMirrorsX2[3] = -2 * infoObject.x1;
      zMirrorsZ[3] = z4;
      zMirrorsP[3].z = z4;
      zMirrorsOP[3] = infoObject.zMirrorZ2OP; // the outer bottom mirror

      break;
    case 1: // 1 = single canonical resonator in x direction
    default:
      xMirrorsN = 2;
      yMirrorsN = 0;
      zMirrorsN = 0;

      for (let i = 0; i < xMirrorsN; i++) {
        xMirrorsY1[i] = y1 + deltaY;
        xMirrorsY2[i] = y2 + deltaY;
        xMirrorsZ1[i] = infoObject.z1;
        xMirrorsZ2[i] = infoObject.z2;
        xMirrorsP[i].y = deltaY;
      }
      xMirrorsX[0] = infoObject.x1;
      xMirrorsP[0].x = infoObject.x1;
      xMirrorsOP[0] = infoObject.xMirrorX1OP;
      xMirrorsX[1] = infoObject.x2;
      xMirrorsP[1].x = infoObject.x2;
      xMirrorsOP[1] = infoObject.xMirrorX2OP;
  }
  infoObject.raytracingSphereShaderMaterial.uniforms.xMirrorsN.value =
    xMirrorsN;
  infoObject.raytracingSphereShaderMaterial.uniforms.yMirrorsN.value =
    yMirrorsN;

  infoObject.raytracingSphereShaderMaterial.uniforms.sphereCentre.value.x =
    infoObject.sphereCentre.x;
  infoObject.raytracingSphereShaderMaterial.uniforms.sphereCentre.value.y =
    infoObject.sphereCentre.y + deltaY;
  infoObject.raytracingSphereShaderMaterial.uniforms.sphereCentre.value.z =
    infoObject.sphereCentre.z;

  // mesh.rotation.y = -Math.atan2(camera.position.z, camera.position.x);
  // mesh.rotation.z = meshRotationZ;

  // raytracingSphereShaderMaterial.uniforms.xMirrorsN.value = xMirrorsN;
  // raytracingSphereShaderMaterial.uniforms.xMirrorsX.value = xMirrorsX;	// {x[0], x[1], ...}; note that we require x[0] <= x[1] <= x[2] ...!
  // raytracingSphereShaderMaterial.uniforms.xMirrorsY1.value = xMirrorsY1;	// {y1[0], y1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.xMirrorsY2.value = xMirrorsY2;	// {y2[0], y2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.xMirrorsZ1.value = xMirrorsZ1;	// {z1[0], z1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.xMirrorsZ2.value = xMirrorsZ2;	// {z2[0], z2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.xMirrorsP.value = xMirrorsP;	// {P[0], P[1], ...}
  // raytracingSphereShaderMaterial.uniforms.xMirrorsOP.value = xMirrorsOP;	// {op[0], op[1], ...}

  // raytracingSphereShaderMaterial.uniforms.yMirrorsN.value = yMirrorsN;
  // raytracingSphereShaderMaterial.uniforms.yMirrorsY.value = yMirrorsY;	// {y[0], y[1], ...}; note that we require y[0] <= y[1] <= y[2] ...!
  // raytracingSphereShaderMaterial.uniforms.yMirrorsX1.value = yMirrorsX1;	// {x1[0], x1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.yMirrorsX2.value = yMirrorsX2;	// {x2[0], x2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.yMirrorsZ1.value = yMirrorsZ1;	// {z1[0], z1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.yMirrorsZ2.value = yMirrorsZ2;	// {z2[0], z2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.yMirrorsP.value = yMirrorsP;	// {P[0], P[1], ...}
  // raytracingSphereShaderMaterial.uniforms.yMirrorsOP.value = yMirrorsOP;	// {op[0], op[1], ...}

  // raytracingSphereShaderMaterial.uniforms.zMirrorsN.value = zMirrorsN;
  // raytracingSphereShaderMaterial.uniforms.zMirrorsZ.value = zMirrorsZ;	// {z[0], z[1], ...}; note that we require z[0] <= z[1] <= z[2] ...!
  // raytracingSphereShaderMaterial.uniforms.zMirrorsX1.value = zMirrorsX1;	// {x1[0], x1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.zMirrorsX2.value = zMirrorsX2;	// {x2[0], x2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.zMirrorsY1.value = zMirrorsY1;	// {y1[0], y1[1], ...}
  // raytracingSphereShaderMaterial.uniforms.zMirrorsY2.value = zMirrorsY2;	// {y2[0], y2[1], ...}
  // raytracingSphereShaderMaterial.uniforms.zMirrorsP.value = zMirrorsP;	// {P[0], P[1], ...}
  // raytracingSphereShaderMaterial.uniforms.zMirrorsOP.value = zMirrorsOP;	// {op[0], op[1], ...}

  infoObject.raytracingSphereShaderMaterial.uniforms.backgroundTexture.value =
    backgroundTexture;

  // create the points on the aperture

  // create basis vectors for the camera's clear aperture
  let viewDirection = new THREE.Vector3();
  let apertureBasisVector1 = new THREE.Vector3();
  let apertureBasisVector2 = new THREE.Vector3();
  infoObject.camera.getWorldDirection(viewDirection);
  viewDirection.normalize();
  // postStatus(`viewDirection.lengthSq() = ${viewDirection.lengthSq()}`);
  // if(counter < 10) console.log(`viewDirection = (${viewDirection.x.toPrecision(2)}, ${viewDirection.y.toPrecision(2)}, ${viewDirection.z.toPrecision(2)})`);

  if (viewDirection.x == 0.0 && viewDirection.y == 0.0) {
    // viewDirection is along z direction
    apertureBasisVector1
      .crossVectors(viewDirection, new THREE.Vector3(1, 0, 0))
      .normalize();
  } else {
    // viewDirection is not along z direction
    apertureBasisVector1
      .crossVectors(viewDirection, new THREE.Vector3(0, 0, 1))
      .normalize();
  }
  apertureBasisVector1
    .crossVectors(THREE.Object3D.DEFAULT_UP, viewDirection)
    .normalize();
  // viewDirection = new THREE.Vector3(0, 0, -1);
  // apertureBasisVector1 = new THREE.Vector3(1, 0, 0);
  apertureBasisVector2
    .crossVectors(viewDirection, apertureBasisVector1)
    .normalize();

  infoObject.raytracingSphereShaderMaterial.uniforms.noOfRays.value =
    infoObject.noOfRays;
  infoObject.raytracingSphereShaderMaterial.uniforms.apertureXHat.value.copy(
    apertureBasisVector1
  );
  infoObject.raytracingSphereShaderMaterial.uniforms.apertureYHat.value.copy(
    apertureBasisVector2
  );
  infoObject.raytracingSphereShaderMaterial.uniforms.viewDirection.value.copy(
    viewDirection
  );
  infoObject.raytracingSphereShaderMaterial.uniforms.apertureRadius.value =
    infoObject.apertureRadius;

  let focusDistance = Math.tan(infoObject.atanFocusDistance);

  if (
    infoObject.raytracingSphereShaderMaterial.uniforms.focusDistance.value !=
    focusDistance
  ) {
    infoObject.raytracingSphereShaderMaterial.uniforms.focusDistance.value =
      focusDistance;
    // GUIParams.'tan<sup>-1</sup>(focus. dist.)'.value = atanFocusDistance;
    focusDistanceControl.setValue(infoObject.atanFocusDistance);
  }

  // (re)create random numbers
  // let i=0;
  // let randomNumbersX = [];
  // let randomNumbersY = [];
  // do {
  // 	// create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
  // 	let x = 2*Math.random()-1;	// random number between -1 and 1
  // 	let y = 2*Math.random()-1;	// random number between -1 and 1
  // 	if(x*x + y*y <= 1) {
  // 		// (x,y) lies within a circle of radius 1
  // 		//  add a new point to the array of points on the aperture
  // 		randomNumbersX.push(apertureRadius*x);
  // 		randomNumbersY.push(apertureRadius*y);
  // 		i++;
  // 	}
  // } while (i < 100);
  // raytracingSphereShaderMaterial.uniforms.randomNumbersX.value = randomNumbersX;
  // raytracingSphereShaderMaterial.uniforms.randomNumbersY.value = randomNumbersY;
}

/** create raytracing phere */
function addRaytracingSphere() {
  // create arrays of random numbers (as GLSL is rubbish at doing random numbers)
  let randomNumbersX = [];
  let randomNumbersY = [];
  // make the first random number 0 in both arrays, meaning the 0th ray starts from the centre of the aperture
  randomNumbersX.push(0);
  randomNumbersY.push(0);
  // fill in the rest of the array with random numbers
  let i = 1;
  do {
    // create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
    let x = 2 * Math.random() - 1; // random number between -1 and 1
    let y = 2 * Math.random() - 1; // random number between -1 and 1
    if (x * x + y * y <= 1) {
      // (x,y) lies within a circle of radius 1
      //  add a new point to the array of points on the aperture
      randomNumbersX.push(x);
      randomNumbersY.push(y);
      i++;
    }
  } while (i < 100);

  // let n2Vectors = [];	// principal points
  // let n2Floats = [];
  // for(i=0; i<mirrorsN2; i++) {
  // 	n2Vectors.push(new THREE.Vector3(0., 0., 0.));
  // 	n2Floats.push(0.);
  // }

  // the sphere surrounding the camera in all directions
  const geometry = new THREE.SphereGeometry(raytracingSphereRadius);
  infoObject.raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
    side: THREE.DoubleSide,
    // wireframe: true,
    uniforms: {
      // the set of mirrors in x planes
      maxTraceLevel: { value: 50 },
      xMirrorsN: { value: 0 },
      xMirrorsX: { value: xMirrorsX }, // {x[0], x[1], ...}; note that we require x[0] <= x[1] <= x[2] ...!
      xMirrorsY1: { value: xMirrorsY1 }, // {y1[0], y1[1], ...}
      xMirrorsY2: { value: xMirrorsY2 }, // {y2[0], y2[1], ...}
      xMirrorsZ1: { value: xMirrorsZ1 }, // {z1[0], z1[1], ...}
      xMirrorsZ2: { value: xMirrorsZ2 }, // {z2[0], z2[1], ...}
      xMirrorsP: { value: xMirrorsP },
      xMirrorsOP: { value: xMirrorsOP },
      // the set of mirrors in y planes
      yMirrorsN: { value: 0 },
      yMirrorsY: { value: yMirrorsY }, // {y[0], y[1], ...}; note that we require y[0] <= y[1] <= y[2] ...!
      yMirrorsX1: { value: yMirrorsX1 }, // {x1[0], x1[1], ...}
      yMirrorsX2: { value: yMirrorsX2 }, // {x2[0], x2[1], ...}
      yMirrorsZ1: { value: yMirrorsZ1 }, // {z1[0], z1[1], ...}
      yMirrorsZ2: { value: yMirrorsZ2 }, // {z2[0], z2[1], ...}
      yMirrorsP: { value: yMirrorsP },
      yMirrorsOP: { value: yMirrorsOP },
      // cylindricalMirrors: { value: true },
      mirrorType: { value: 1 },
      reflectionCoefficient: { value: 0.9 },
      sphereCentre: { value: new THREE.Vector3(0, 0, 0) },
      sphereRadius: { value: sphereRadius },
      sphereHeight: { value: sphereHeight },
      showSphere: { value: true },
      outerRadius: { value: outerRadius },
      outerHeightNegative: { value: outerHeightNegative },
      outerHeightPositive: { value: outerHeightPositive },
      outerPhaseShift: { value: outerPhaseShift },
      innerRadius: { value: innerRadius },
      innerHeightNegative: { value: innerHeightNegative },
      innerHeightPositive: { value: innerHeightPositive },
      innerPhaseShift: { value: innerPhaseShift },
      showInnerCylinder: { value: true },
      showOuterCylinder: { value: true },
      showCloak: { value: true },
      backgroundTexture: { value: backgroundTexture },
      focusDistance: { value: 10.0 },
      apertureXHat: { value: new THREE.Vector3(1, 0, 0) },
      apertureYHat: { value: new THREE.Vector3(0, 1, 0) },
      apertureRadius: { value: infoObject.apertureRadius },
      randomNumbersX: { value: randomNumbersX },
      randomNumbersY: { value: randomNumbersY },
      noOfRays: { value: 1 },
      viewDirection: { value: new THREE.Vector3(0, 0, -1) },
      keepVideoFeedForward: { value: true },
    },
    vertexShader: `varying vec3 intersectionPoint;

    void main()	{
      // projectionMatrix, modelViewMatrix, position -> passed in from Three.js
      intersectionPoint = (modelMatrix * vec4(position, 1.0)).xyz;	// position.xyz;
      
        gl_Position = projectionMatrix
        * modelViewMatrix
        * vec4(position, 1.0);
    }`,
    fragmentShader: `precision highp float;

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
      float r,	// sphere radius
      out vec3 intersectionPosition,
      out float intersectionDistance
    ) {
      c = c + vec3(0, y , 0);
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
      if( showCloak && showInnerCylinder && findNearestIntersectionWithCylinder(s, d, sphereCentre, innerRadius, innerHeightNegative, innerHeightPositive, startIntersectionObjectIndex == 0, intersectionPosition, intersectionDistance, intersectionNormal) ) {
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
      if( showCloak && showOuterCylinder && findNearestIntersectionWithCylinder(s, d, sphereCentre, outerRadius, outerHeightNegative, outerHeightPositive, startIntersectionObjectIndex == 1, intersectionPosition, intersectionDistance, intersectionNormal) ) {
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
      if( showSphere && findNearestIntersectionWithSphere(s, d, sphereCentre, sphereHeight, sphereRadius, intersectionPosition, intersectionDistance) ) {
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
    }`
  });

  raytracingSphere = new THREE.Mesh(
    geometry,
    infoObject.raytracingSphereShaderMaterial
  );
  scene.add(raytracingSphere);
}

// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_animation_skinning_additive_blending.html
// the best thing to do would be to move the whole createGUI() into its own file because this is a mess
// todo: add presets for cloak settings
function createGUI() {
  // const
  gui = new GUI();
  // gui.hide();
  // GUIMesh = new HTMLMesh( gui.domElement );	// placeholder

  GUIParams = {
    noOfReflections:
      infoObject.raytracingSphereShaderMaterial.uniforms.maxTraceLevel.value -
      2,
    "Horiz. FOV (&deg;)": infoObject.fovScreen,
    "Aperture radius": infoObject.apertureRadius,
    "tan<sup>-1</sup>(focus. dist.)": infoObject.atanFocusDistance,
    "No of rays": infoObject.noOfRays,
    autofocus: function () {
      autofocus = !autofocus;
      autofocusControl.name("Autofocus: " + (autofocus ? "On" : "Off"));
      focusDistanceControl.disable(autofocus);
    }, // (autofocus?'On':'Off'),
    // 'Autofocus': autofocus,
    "Point forward (in -<b>z</b> direction)": () =>
      pointForward(infoObject, orbitControls),
    "Show/hide info": () => {
      toggleInfoVisibility(infoObject);
    },
    vrControlsVisible: function () {
      GUIMesh.visible = !GUIMesh.visible;
      vrControlsVisibleControl.name(guiMeshVisible2String(GUIMesh));
    },
    background: () => {
      infoObject.background = (infoObject.background + 1) % 5;
      backgroundTexture = loadBackgroundImage(
        infoObject.background,
        backgroundTexture
      );
      backgroundControl.name(background2String(infoObject.background));
    },
    sphereRadius: sphereRadius,
    sphereHeight: sphereHeight,
    outerRadius: outerRadius,
    outerHeightNegative: outerHeightNegative,
    outerHeightPositive: outerHeightPositive,
    outerPhaseShift: outerPhaseShift,
    innerRadius: innerRadius,
    innerHeightNegative: innerHeightNegative,
    innerHeightPositive: innerHeightPositive,
    innerPhaseShift: innerPhaseShift,
    sphereCentreX: infoObject.sphereCentre.x,
    sphereCentreY: infoObject.sphereCentre.y,
    sphereCentreZ: infoObject.sphereCentre.z,
    showSphere: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.showSphere.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.showSphere.value;
      showSphereControl.name(
        showSphere2String(infoObject.raytracingSphereShaderMaterial)
      );
    },
    showCloak: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.showCloak.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.showCloak.value;
      showCloakControl.name(
        showCloak2String(infoObject.raytracingSphereShaderMaterial)
      );
    },

    showInnerCylinder: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.showInnerCylinder.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.showInnerCylinder
          .value;
      showInnerCylinderControl.name(
        showInnerCylinder2String(infoObject.raytracingSphereShaderMaterial)
      );
    },
    showOuterCylinder: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.showOuterCylinder.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.showOuterCylinder
          .value;
      showOuterCylinderControl.name(
        showOuterCylinder2String(infoObject.raytracingSphereShaderMaterial)
      );
    },
    resonatorType: function () {
      infoObject.resonatorType = (infoObject.resonatorType + 1) % 4;
      resonatorTypeControl.name(resonatorType2String(infoObject.resonatorType));
      enableDisableResonatorControls();
      if (infoObject.resonatorType == 3) {
        // Penrose cavity
        infoObject.zMirrorZ1OP = Math.max(1, infoObject.zMirrorZ1OP);
        infoObject.zMirrorZ2OP = Math.max(1, infoObject.zMirrorZ2OP);
        opz0Control.setValue(infoObject.zMirrorZ1OP);
        opz1Control.setValue(infoObject.zMirrorZ2OP);
      }
      // createGUI();
      // opz0Control.disable( resonatorType == 0 );
      // opz1Control.disable( resonatorType == 0 );
      // z0Control.disable( resonatorType == 0 );
      // z1Control.disable( resonatorType == 0 );
    },
    mirrorType: function () {
      infoObject.raytracingSphereShaderMaterial.uniforms.mirrorType.value =
        (infoObject.raytracingSphereShaderMaterial.uniforms.mirrorType.value +
          1) %
        2;
      mirrorTypeControl.name(
        mirrorType2String(
          infoObject.raytracingSphereShaderMaterial.uniforms.mirrorType.value
        )
      );
    },
    // optical powers
    // opx1: infoObject.xMirrorX1OP,
    // opx2: infoObject.xMirrorX2OP,
    // opz1: infoObject.zMirrorZ1OP,
    // opz2: infoObject.zMirrorZ2OP,
    x1: infoObject.x1,
    // x2: infoObject.x2,
    // z1: infoObject.z1,
    // z2: infoObject.z2,
    resonatorY: infoObject.resonatorY,
    cylindricalMirrors: function () {
      infoObject.raytracingSphereShaderMaterial.uniforms.cylindricalMirrors.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.cylindricalMirrors
          .value;
      cylindricalMirrorsControl.name(cylindricalMirrors2String());
    },
    // reflectionCoefficient9s: -Math.log10(1-raytracingSphereShaderMaterial.uniforms.reflectionCoefficient.value),
    reflectionLossDB:
      10 *
      Math.log10(
        1 -
          infoObject.raytracingSphereShaderMaterial.uniforms
            .reflectionCoefficient.value
      ),
    makeEyeLevel: function () {
      infoObject.resonatorY = infoObject.camera.position.y;
      resonatorYControl.setValue(infoObject.resonatorY);
    },
    // meshRotX: meshRotationX,
    // meshRotY: meshRotationY,
    // meshRotZ: meshRotationZ
  };

  gui
    .add(GUIParams, "noOfReflections", 0, 200, 1)
    .name("Max. reflections")
    .onChange((r) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.maxTraceLevel.value =
        r + 2;
    });
  gui
    .add(GUIParams, "outerPhaseShift", 0, 1, 0.1)
    .name("Outer hologram phase shift")
    .onChange((outer_Shift) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerPhaseShift.value =
        outer_Shift;
      console.log(outer_Shift);
    });
  gui
    .add(GUIParams, "innerPhaseShift", -1, 0, 0.1)
    .name("Inner hologram phase shift")
    .onChange((inner_Shift) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerPhaseShift.value =
        inner_Shift;
    });
  // mirrorTypeControl = gui
  //   .add(GUIParams, "mirrorType")
  //   .name(
  //     mirrorType2String(
  //       infoObject.raytracingSphereShaderMaterial.uniforms.mirrorType.value
  //     )
  //   );
  // cylindricalMirrorsControl = gui.add( GUIParams, 'cylindricalMirrors' ).name( cylindricalMirrors2String() );
  // gui
  //   .add(GUIParams, "x1", -10, -0.1, 0.001)
  //   .name("<i>x</i><sub>1</sub>")
  //   .onChange((x) => {
  //     infoObject.x1 = x;
  //   });

  resonatorYControl = gui
    .add(GUIParams, "resonatorY", 0, 3, 0.001)
    .name("<i>y</i><sub>resonator</sub>")
    .onChange((y) => {
      infoObject.resonatorY = y;
      refreshInfo(infoObject);
    });

  // gui.add(GUIParams, "makeEyeLevel").name("Move resonator to eye level");
  // // gui.add( GUIParams, 'reflectionCoefficient9s', 0, 3, 0.1 ).name( '<div class="tooltip">Nines(<i>R</i>)<span class="tooltiptext">The number of <a href="https://en.m.wikipedia.org/wiki/Nines_(notation)">nines</a><br>in the reflection<br>coefficient, <i>R</i>.<br>E.g. Nines(0.99) = 2.</span></div> ' ).onChange( (l) => { raytracingSphereShaderMaterial.uniforms.reflectionCoefficient.value = 1-Math.pow(10, -l); } );
  // gui
  //   .add(GUIParams, "reflectionLossDB", -30, 0, 0.1)
  //   .name("Refl. loss (dB)")
  //   .onChange((l) => {
  //     infoObject.raytracingSphereShaderMaterial.uniforms.reflectionCoefficient.value =
  //       1 - Math.pow(10, 0.1 * l);
  //   });
  // remove these for the moment
  // gui.add( GUIParams, 'sphereCentreX', -5, 5 ).name( "<i>x</i><sub>sphere</sub>" ).onChange( (x) => { sphereCentre.x = x; } );
  // gui.add( GUIParams, 'sphereCentreY',  0, 5 ).name( "<i>y</i><sub>sphere</sub>" ).onChange( (y) => { sphereCentre.y = y; } );
  // gui.add( GUIParams, 'sphereCentreZ', -5, 5 ).name( "<i>z</i><sub>sphere</sub>" ).onChange( (z) => { sphereCentre.z = z; } );

  const sphereFolder = gui.addFolder("Sphere Controls");

  showSphereControl = sphereFolder
    .add(GUIParams, "showSphere")
    .name(showSphere2String(infoObject.raytracingSphereShaderMaterial));

  sphereFolder
    .add(GUIParams, "sphereRadius", 0, 1)
    .name("<i>r</i><sub>sphere</sub>")
    .onChange((r) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.sphereRadius.value = r;
    });
  sphereFolder
    .add(GUIParams, "sphereHeight", -1, 1, 0.05)
    .name("<i>h</i><sub>sphere</sub>")
    .onChange((h_sphere) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.sphereHeight.value =
        h_sphere;
    });

  const cloakFolder = gui.addFolder("Axicon Cloak Controls");

  showCloakControl = cloakFolder
    .add(GUIParams, "showCloak")
    .name(showCloak2String(infoObject.raytracingSphereShaderMaterial));

  showOuterCylinderControl = cloakFolder
    .add(GUIParams, "showOuterCylinder")
    .name(showOuterCylinder2String(infoObject.raytracingSphereShaderMaterial));

  cloakFolder
    .add(GUIParams, "innerRadius", 0, 1)
    .name("<i>r</i><sub>innner</sub>")
    .onChange((r_inner) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerRadius.value =
        r_inner;
    });
  cloakFolder
    .add(GUIParams, "outerRadius", 0, 1)
    .name("<i>r</i><sub>outer</sub>")
    .onChange((r_outer) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerRadius.value =
        r_outer;
    });

  cloakFolder
    .add(GUIParams, "outerHeightNegative", -1, 0, 0.1)
    .name("<i>h</i><sub>-outer</sub>")
    .onChange((h_outer_neg) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerHeightNegative.value =
        h_outer_neg;
    });
  cloakFolder
    .add(GUIParams, "outerHeightPositive", 0, 1, 0.1)
    .name("<i>h</i><sub>+outer</sub>")
    .onChange((h_outer_pos) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerHeightPositive.value =
        h_outer_pos;
    });

  showInnerCylinderControl = cloakFolder
    .add(GUIParams, "showInnerCylinder")
    .name(showInnerCylinder2String(infoObject.raytracingSphereShaderMaterial));
  cloakFolder
    .add(GUIParams, "innerRadius", 0, 1)
    .name("<i>r</i><sub>innner</sub>")
    .onChange((r_inner) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerRadius.value =
        r_inner;
    });

  cloakFolder
    .add(GUIParams, "innerHeightNegative", -1, 0, 0.1)
    .name("<i>h</i><sub>-inner</sub>")
    .onChange((h_inner_neg) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerHeightNegative.value =
        h_inner_neg;
    });
  cloakFolder
    .add(GUIParams, "innerHeightPositive", 0, 1, 0.1)
    .name("<i>h</i><sub>+inner</sub>")
    .onChange((h_inner_pos) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerHeightPositive.value =
        h_inner_pos;
    });

  gui.add(GUIParams, "Point forward (in -<b>z</b> direction)");
  backgroundControl = gui
    .add(GUIParams, "background")
    .name(background2String(infoObject.background));
  // gui.add( GUIParams, 'meshRotX', -Math.PI, Math.PI ).name('Rot x').onChange( (a) => { meshRotationX = a; })
  // gui.add( GUIParams, 'meshRotY', -Math.PI, Math.PI ).name('Rot y').onChange( (a) => { meshRotationY = a; })
  // gui.add( GUIParams, 'meshRotZ', -Math.PI, Math.PI ).name('Rot z').onChange( (a) => { meshRotationZ = a; })

  // const folderVirtualCamera = gui.addFolder( 'Virtual camera' );
  gui.add(GUIParams, "Horiz. FOV (&deg;)", 1, 170, 1).onChange((fov) => {
    screenChanged(renderer, infoObject.camera, fov);
    infoObject.fovScreen = fov;
  });
  gui.add(GUIParams, "Aperture radius", 0.0, 1.0, 0.01).onChange((r) => {
    infoObject.apertureRadius = r;
  });
  // autofocusControl = gui.add( GUIParams, 'autofocus' ).name( 'Autofocus: ' + (autofocus?'On':'Off') );
  // gui.add( GUIParams, 'Autofocus' ).onChange( (b) => { autofocus = b; focusDistanceControl.disable(autofocus); } );
  focusDistanceControl = gui
    .add(
      GUIParams,
      "tan<sup>-1</sup>(focus. dist.)",
      //Math.atan(0.1),
      0.01, // -0.5*Math.PI,	// allow only positive focussing distances
      0.5 * Math.PI,
      0.0001
    )
    .onChange((a) => {
      infoObject.atanFocusDistance = a;
    });
  focusDistanceControl.disable(autofocus);
  // focusDistanceControl = gui.add( GUIParams, 'tan<sup>-1</sup>(focus. dist.)',
  // 	//Math.atan(0.1),
  // 	-0.5*Math.PI,
  // 	0.5*Math.PI,
  // 	0.001
  // ).onChange( (a) => { atanFocusDistance = a; } );
  // folderVirtualCamera.add( atanFocusDistance, 'atan focus dist', -0.5*Math.PI, +0.5*Math.PI ).listen();
  gui.add(GUIParams, "No of rays", 1, 100, 1).onChange((n) => {
    infoObject.noOfRays = n;
  });

  if (renderer.xr.enabled) {
    vrControlsVisibleControl = gui.add(GUIParams, "vrControlsVisible");
  }
  // folderVirtualCamera.close();

  // const folderSettings = gui.addFolder( 'Other controls' );
  // // folderSettings.add( params, 'Video feed forward' ).onChange( (b) => { raytracingSphereShaderMaterial.uniforms.keepVideoFeedForward.value = b; } );
  // // folderSettings.add( params, 'Lenslet type', { 'Ideal thin': true, 'Phase hologram': false } ).onChange( (t) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = t; });
  // // folderSettings.add( params, 'Ideal lenses').onChange( (b) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = b; } );
  // folderSettings.add( params, 'Show/hide info');
  // folderSettings.close();

  // enableDisableResonatorControls();

  // create the GUI mesh at the end to make sure that it includes all controls
  GUIMesh = new HTMLMesh(gui.domElement);
  GUIMesh.visible = false;
  vrControlsVisibleControl.name(guiMeshVisible2String(GUIMesh)); // this can be called only after GUIMesh has been created
}

//at this point does nothing?

// function cylindricalMirrors2String() {
// 	return (raytracingSphereShaderMaterial.uniforms.cylindricalMirrors.value?'Cylindrical mirrors':'Spherical mirrors');
// 	// return (raytracingSphereShaderMaterial.uniforms.cylindricalMirrors.value?'<div class="tooltip">Cylindrical mirrors<span class="tooltiptext">Cylindrical mirrors,<br>simulated as planar,<br>reflective, ideal<br>thin cylindrical<br>lenses</span></div>':'<div class="tooltip">Spherical mirrors<span class="tooltiptext">Spherical mirrors,<br>simulated as planar,<br>reflective, ideal<br>thin lenses</span></div>');
// }
//Do I move this function to the InfoFunctions.js or is this useless?

function initMirrors() {
  xMirrorsN = 2;
  // initialise all the elements to default values
  xMirrorsX = []; // {x[0], x[1], ...}; note that we require x[0] <= x[1] <= x[2] ...!
  xMirrorsY1 = []; // {y1[0], y1[1], ...}
  xMirrorsY2 = []; // {y2[0], y2[1], ...}
  xMirrorsZ1 = []; // {z1[0], z1[1], ...}
  xMirrorsZ2 = []; // {z2[0], z2[1], ...}
  xMirrorsP = []; // principal points
  xMirrorsOP = []; // optical powers
  for (let i = 0; i < mirrorsN2; i++) {
    xMirrorsX.push(0.0);
    xMirrorsY1.push(y1 + infoObject.resonatorY);
    xMirrorsY2.push(y2 + infoObject.resonatorY);
    xMirrorsZ1.push(infoObject.z1);
    xMirrorsZ2.push(infoObject.z2);
    xMirrorsP.push(new THREE.Vector3(0, infoObject.resonatorY, 0));
    xMirrorsOP.push(-1);
  }
  // set the actual values where those differ from the default ones
  xMirrorsX[0] = infoObject.x1;
  xMirrorsP[0].x = infoObject.x1;
  // xMirrorsX[1] = 0.2; xMirrorsP[1].x = 0.2; xMirrorsY2[1] = 2; xMirrorsZ2[1] = 0;
  xMirrorsX[1] = infoObject.x2;
  xMirrorsP[1].x = infoObject.x2;

  yMirrorsN = 0;
  // initialise all the elements to default values
  yMirrorsY = []; // {y[0], y[1], ...}; note that we require y[0] <= y[1] <= y[2] ...!
  yMirrorsX1 = []; // {x1[0], x1[1], ...}
  yMirrorsX2 = []; // {x2[0], x2[1], ...}
  yMirrorsZ1 = []; // {z1[0], z1[1], ...}
  yMirrorsZ2 = []; // {z2[0], z2[1], ...}
  yMirrorsP = []; // principal points
  yMirrorsOP = []; // optical powers
  for (let i = 0; i < mirrorsN2; i++) {
    yMirrorsY.push(0.0);
    yMirrorsX1.push(infoObject.x1);
    yMirrorsX2.push(infoObject.x2);
    yMirrorsZ1.push(infoObject.z1);
    yMirrorsZ2.push(infoObject.z2);
    yMirrorsP.push(new THREE.Vector3(0, infoObject.resonatorY, 0));
    yMirrorsOP.push(0);
  }
  // set the actual values where those differ from the default ones
  yMirrorsY[0] = y1 + infoObject.resonatorY;
  yMirrorsP[0].y = y1 + infoObject.resonatorY;
  yMirrorsY[1] = y2 + infoObject.resonatorY;
  yMirrorsP[1].y = y2 + infoObject.resonatorY;

  zMirrorsN = 2;
  // initialise all the elements to default values
  zMirrorsZ = []; // {z[0], z[1], ...}; note that we require z[0] <= z[1] <= z[2] ...!
  zMirrorsX1 = []; // {x1[0], x1[1], ...}
  zMirrorsX2 = []; // {x2[0], x2[1], ...}
  zMirrorsY1 = []; // {y1[0], y1[1], ...}
  zMirrorsY2 = []; // {y2[0], y2[1], ...}
  zMirrorsP = []; // principal points
  zMirrorsOP = []; // optical powers
  for (let i = 0; i < mirrorsN2; i++) {
    zMirrorsZ.push(0.0);
    zMirrorsX1.push(infoObject.x1);
    zMirrorsX2.push(infoObject.x2);
    zMirrorsY1.push(y1 + infoObject.resonatorY);
    zMirrorsY2.push(y2 + infoObject.resonatorY);
    zMirrorsP.push(new THREE.Vector3(0, infoObject.resonatorY, 0));
    zMirrorsOP.push(0.1);
  }
  // set the actual values where those differ from the default ones
  zMirrorsZ[0] = infoObject.z1;
  zMirrorsP[0].z = infoObject.z1;
  zMirrorsZ[1] = infoObject.z2;
  zMirrorsP[1].z = infoObject.z2;
}

function addXRInteractivity() {
  // see https://github.com/mrdoob/three.js/blob/master/examples/webxr_vr_sandbox.html

  // the two hand controllers

  const geometry = new THREE.BufferGeometry();
  geometry.setFromPoints([
    new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, 0, -5),
  ]);

  const controller1 = renderer.xr.getController(0);
  controller1.add(new THREE.Line(geometry));
  scene.add(controller1);

  const controller2 = renderer.xr.getController(1);
  controller2.add(new THREE.Line(geometry));
  scene.add(controller2);

  //

  const controllerModelFactory = new XRControllerModelFactory();

  const controllerGrip1 = renderer.xr.getControllerGrip(0);
  controllerGrip1.add(
    controllerModelFactory.createControllerModel(controllerGrip1)
  );
  scene.add(controllerGrip1);

  const controllerGrip2 = renderer.xr.getControllerGrip(1);
  controllerGrip2.add(
    controllerModelFactory.createControllerModel(controllerGrip2)
  );
  scene.add(controllerGrip2);

  //

  const group = new InteractiveGroup(renderer, infoObject.camera);
  group.listenToPointerEvents(renderer, infoObject.camera);
  group.listenToXRControllerEvents(controller1);
  group.listenToXRControllerEvents(controller2);
  scene.add(group);

  // place this below the resonator
  // GUIMesh = new HTMLMesh( gui.domElement );
  GUIMesh.position.x = 0;
  GUIMesh.position.y = infoObject.resonatorY - 1.5;
  GUIMesh.position.z = -0.4;
  GUIMesh.rotation.x = -Math.PI / 4;
  GUIMesh.scale.setScalar(2);
  group.add(GUIMesh);
}

function createVideoFeeds() {
  // create the video stream for the user-facing camera first, as some devices (such as my iPad), which have both cameras,
  // but can (for whatever reason) only have a video feed from one at a time, seem to go with the video stream that was
  // created last, and as the standard view is looking "forward" it is preferable to see the environment-facing camera.
  videoFeedU = document.getElementById("videoFeedU");

  // see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
  if (navigator.mediaDevices && navigator.mediaDevices.getUserMedia) {
    // user-facing camera
    const constraintsU = {
      video: {
        // 'deviceId': cameraId,	// this could be the device ID selected
        width: { ideal: 1280 }, // {ideal: 10000},
        // height: {ideal: 10000},
        facingMode: { ideal: "user" },
        // aspectRatio: { exact: width / height }
      },
    };
    navigator.mediaDevices
      .getUserMedia(constraintsU)
      .then(function (stream) {
        // apply the stream to the video element used in the texture
        videoFeedU.srcObject = stream;
        videoFeedU.play();

        videoFeedU.addEventListener("playing", () => {
          aspectRatioVideoFeedU =
            videoFeedU.videoWidth / videoFeedU.videoHeight;
          updateUniforms();
          postStatus(
            `User-facing(?) camera resolution ${videoFeedU.videoWidth} &times; ${videoFeedU.videoHeight}`
          );
        });
      })
      .catch(function (error) {
        postStatus(
          `Unable to access user-facing camera/webcam (Error: ${error})`
        );
      });
  } else {
    postStatus(
      "MediaDevices interface, which is required for video streams from device cameras, not available."
    );
  }

  videoFeedE = document.getElementById("videoFeedE");

  // see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
  if (navigator.mediaDevices && navigator.mediaDevices.getUserMedia) {
    // environment-facing camera
    const constraintsE = {
      video: {
        // 'deviceId': cameraId,	// this could be the device ID selected
        width: { ideal: 1280 }, // {ideal: 10000},
        // height: {ideal: 10000},
        facingMode: { ideal: "environment" },
        // aspectRatio: { exact: width / height }
      },
    };
    navigator.mediaDevices
      .getUserMedia(constraintsE)
      .then(function (stream) {
        // apply the stream to the video element used in the texture
        videoFeedE.srcObject = stream;
        videoFeedE.play();

        videoFeedE.addEventListener("playing", () => {
          aspectRatioVideoFeedE =
            videoFeedE.videoWidth / videoFeedE.videoHeight;
          updateUniforms();
          postStatus(
            `Environment-facing(?) camera resolution ${videoFeedE.videoWidth} &times; ${videoFeedE.videoHeight}`
          );
        });
      })
      .catch(function (error) {
        postStatus(
          `Unable to access environment-facing camera/webcam (Error: ${error})`
        );
      });
  } else {
    postStatus(
      "MediaDevices interface, which is required for video streams from device cameras, not available."
    );
  }
}

function addEventListenersEtc() {
  // handle device orientation
  // window.addEventListener("deviceorientation", handleOrientation, true);

  // handle window resize
  window.addEventListener(
    "resize",
    () => {
      onWindowResize(renderer, infoObject.camera, infoObject.fovScreen);
    },
    false
  );

  // handle screen-orientation (landscape/portrait) change
  screen.orientation.addEventListener("change", recreateVideoFeeds);

  // share button functionality
  document.getElementById("takePhotoButton").addEventListener("click", () => {
    takePhoto(storedPhoto, renderer, infoObject);
    //storedPhoto = takePhoto(storedPhoto, renderer, infoObject);
  });

  // toggle fullscreen button functionality
  document
    .getElementById("fullscreenButton")
    .addEventListener("click", toggleFullscreen);

  // info button functionality
  document.getElementById("infoButton").addEventListener("click", () => {
    toggleInfoVisibility(infoObject);
  });

  // back button functionality
  document
    .getElementById("backButton")
    .addEventListener("click", () => showLivePhoto(renderer, gui, infoObject));
  document.getElementById("backButton").style.visibility = "hidden";

  // share button
  document.getElementById("shareButton").addEventListener("click", share);
  document.getElementById("shareButton").style.visibility = "hidden";
  if (!navigator.share)
    document.getElementById("shareButton").src = "./shareButtonUnavailable.png";
  // if(!(navigator.share)) document.getElementById('shareButton').style.opacity = 0.3;

  // delete button
  document.getElementById("deleteButton").addEventListener("click", () => {
    deleteStoredPhoto(storedPhoto, renderer, gui, infoObject);
  });
  document.getElementById("deleteButton").style.visibility = "hidden";

  // hide the thumbnail for the moment
  document
    .getElementById("storedPhotoThumbnail")
    .addEventListener("click", () =>
      showStoredPhoto(renderer, gui, infoObject)
    );

  document.getElementById("storedPhotoThumbnail").style.visibility = "hidden";

  document
    .getElementById("storedPhoto")
    .addEventListener("click", () => showLivePhoto(renderer, gui, infoObject));
  document.getElementById("storedPhoto").style.visibility = "hidden";
  // showingStoredPhoto = false;
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeeds() {
  // stop current video streams...
  videoFeedE.srcObject.getTracks().forEach(function (track) {
    track.stop();
  });
  videoFeedU.srcObject.getTracks().forEach(function (track) {
    track.stop();
  });

  // ... and re-create new ones, hopefully of the appropriate size
  createVideoFeeds();
}

function addDragControls() {
  let objects = [];
  objects.push(GUIMesh);

  dragControls = new DragControls(
    objects,
    infoObject.camera,
    renderer.domElement
  );

  // add event listener to highlight dragged objects
  dragControls.addEventListener("dragstart", function (event) {
    event.object.material.emissive.set(0xaaaaaa);
  });

  dragControls.addEventListener("dragend", function (event) {
    event.object.material.emissive.set(0x000000);
  });
}

async function share() {
  try {
    fetch(storedPhoto)
      .then((response) => response.blob())
      .then((blob) => {
        const file = new File(
          [blob],
          infoObject.storedPhotoDescription + ".png",
          {
            type: blob.type,
          }
        );

        // Use the Web Share API to share the screenshot
        if (navigator.share) {
          navigator.share({
            title: infoObject.storedPhotoDescription,
            text: infoObject.storedPhotoInfoString,
            files: [file],
          });
        } else {
          postStatus("Sharing is not supported by this browser.");
        }
      })
      .catch((error) => {
        console.error("Error:", error);
        postStatus(`Error: ${error}`);
      });
  } catch (error) {
    console.error("Error:", error);
  }
}
