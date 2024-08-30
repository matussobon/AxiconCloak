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

import { addOrbitControls, pointForward } from "./camerafunctions.js";
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
  showLens2String,
  showInnerCylinder2String,
  showOuterCylinder2String,
} from "./infofunctions.js";
import {
  postStatus,
  createStatus,
  screenChanged,
  onWindowResize,
} from "./statusfunctions.js";
import { loadBackgroundImage, toggleFullscreen } from "./cosmeticsfunctions.js";
import {
  takePhoto,
  showStoredPhoto,
  showLivePhoto,
  deleteStoredPhoto,
} from "./photofunctions.js";


const NONE = -1;

// surface types
const SURFACE_TYPE_LENS = 0;
const SURFACE_TYPE_COLOR = 1;

// lens types
const LENS_TYPE_IDEAL = 0;
const LENS_TYPE_HOLOGRAM = 1;

// this works fine when running the app locally
// import vertexShaderCode from "./vertex_shader_test.glsl";
// import fragmentShaderCode from "./fragment_shader_test3.glsl";
// a vite.config.mjs has to be included
//https://www.youtube.com/watch?v=RDughHM9qoE
//https://www.npmjs.com/package/vite-plugin-glsl

// random reddit comment came in clutch
// this works fine both locally and when deployed on github
const responseFragment = await fetch("./fragment_shader.glsl");
const fragmentShaderCode = await responseFragment.text();

const responseVertex = await fetch("./vertex_shader.glsl");
const vertexShaderCode = await responseVertex.text();

let scene;
let renderer;
let backgroundTexture;
let orbitControls;
let dragControls;
let raytracingSphere;

let yShift = 0;

let y1 = -0.5;
let y2 = 0.5;

let sphereRadius = 0.1;
let sphereHeight = 0;

let outerRadius = 0.3;
let outerHeightNegative = -0.1;
let outerHeightPositive = 0.1;
let outerYcoord = 0;

let innerRadius = 0.1;
let innerHeightNegative = -0.2;
let innerHeightPositive = 0.2;
let innerYcoord = 0;

let phaseShift = 0.6;
let rotAngle = 0;

let raytracingSphereRadius = 100.0;

let autofocus = false;

// the menu
let gui;
let GUIParams;
let autofocusControl,
  focusDistanceControl,
  resonatorYControl,
  cylindricalMirrorsControl,
  backgroundControl,
  vrControlsVisibleControl,
  showSphereControl;

let showCloakControl;
let showInnerCylinderControl;
let showOuterCylinderControl;

let showLensControl;

let GUIMesh;
let showGUIMesh;
// let meshRotationX = -Math.PI/4, meshRotationY = 0, meshRotationZ = 0;

let lensSurface0, lensSurface1;
let lensSurfaces = [];

let rectangle0, rectangle1;
let rectangles = [];

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

  // rectangle0 = {
  //   corner: new THREE.Vector3(-0.5, 0, 0),
  //   uSpanVector: new THREE.Vector3(1, 0, 0),
  //   vSpanVector: new THREE.Vector3(0, 1, 0),
  //   uSize: 1.0,
  //   vSize: 1.0,
  //   surfaceType: SURFACE_TYPE_LENS,
  //   surfaceIndex: 0,
  // };

  // rectangle1 = {
  //   corner: new THREE.Vector3(-0.5, 0, 0),
  //   uSpanVector: new THREE.Vector3(1, 0, 0),
  //   vSpanVector: new THREE.Vector3(0, -1, 0),
  //   uSize: 1,
  //   vSize: 1,
  //   surfaceType: SURFACE_TYPE_LENS,
  //   surfaceIndex: 1,
  // };

  // rectangles.push(rectangle0, rectangle1);

  // lensSurface0 = {
  //   principalPoint: new THREE.Vector3(0, 0.5, 0),
  //   opticalAxisDirection: new THREE.Vector3(0, 0, 1),
  //   focalLength: 10,
  //   transmissionCoefficient: 0.95,
  //   lensType: LENS_TYPE_IDEAL
  // };

  // lensSurface1 = {
  //   principalPoint: new THREE.Vector3(0, -0.5, 0),
  //   opticalAxisDirection: new THREE.Vector3(0, 0, 1),
  //   focalLength: 1,
  //   transmissionCoefficient: 0.95,
  //   lensType: LENS_TYPE_IDEAL
  // };

  // lensSurfaces.push(lensSurface0, lensSurface1);

  // addLensFan();

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

function addLensFan() {
  let i;
  let vis = [true, true, true];
  for(i=0; i<3; i++) {
    let corner = new THREE.Vector3(0, 0, 0);
    let u = new THREE.Vector3(0, 1, 0);
    let v = new THREE.Vector3(Math.sin(1*i), 0, Math.cos(1*i));
    let p = (new THREE.Vector3(0, 0, 0)).copy(corner).addScaledVector(u, 0.5).addScaledVector(v, 0.5);
    let a = (new THREE.Vector3(0, 0, 0)).crossVectors(u, v);
    let rectangleTemp = {
      visible: vis[i],
      corner: corner,
      uSpanVector: u,
      vSpanVector: v,
      uSize: 1,
      vSize: 1,
      surfaceType: SURFACE_TYPE_LENS,
      surfaceIndex: lensSurfaces.length,
    };
    rectangles.push(rectangleTemp);

    let lensSurfaceTemp = {
      principalPoint: p,
      opticalAxisDirection: a,
      focalLength: i+1,
      transmissionCoefficient: 0.95,
      lensType: LENS_TYPE_IDEAL
    };
    lensSurfaces.push(lensSurfaceTemp);
  }
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

  GUIMesh.position.y = deltaY - 1;
  infoObject.raytracingSphereShaderMaterial.uniforms.backgroundTexture.value =
    backgroundTexture;

  // create the points on the aperture

  // create basis vectors for the camera's clear aperture
  let viewDirection = new THREE.Vector3();
  let apertureBasisVector1 = new THREE.Vector3();
  let apertureBasisVector2 = new THREE.Vector3();
  infoObject.camera.getWorldDirection(viewDirection);
  viewDirection.normalize();

  // this bit should be useless?
  // if (viewDirection.x == 0.0 && viewDirection.y == 0.0) {
  //   // viewDirection is along z direction
  //   apertureBasisVector1
  //     .crossVectors(viewDirection, new THREE.Vector3(1, 0, 0))
  //     .normalize();
  // } else {
  //   // viewDirection is not along z direction
  //   apertureBasisVector1
  //     .crossVectors(viewDirection, new THREE.Vector3(0, 0, 1))
  //     .normalize();
  // }

  apertureBasisVector1
    .crossVectors(THREE.Object3D.DEFAULT_UP, viewDirection)
    .normalize();

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
  }

  rectangles = [];
  lensSurfaces = [];
  addLensFan();
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

  let colorRed = {
    color: new THREE.Vector4(1, 0, 0, 1),
  };

  let colorGreen = {
    color: new THREE.Vector4(0, 1, 0, 1),
  };

  rectangles = [];
  lensSurfaces = [];
  addLensFan();

  // the sphere surrounding the camera in all directions
  const geometry = new THREE.SphereGeometry(raytracingSphereRadius);
  infoObject.raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
    side: THREE.DoubleSide,
    // wireframe: true,
    uniforms: {
      // the set of mirrors in x planes
      maxTraceLevel: { value: 50 },
      // cylindricalMirrors: { value: true },
      mirrorType: { value: 1 },
      reflectionCoefficient: { value: 0.9 },
      sphereCentre: { value: new THREE.Vector3(0, 0, 0) },
      sphereRadius: { value: sphereRadius },
      sphereHeight: { value: sphereHeight },
      showSphere: { value: false },
      outerRadius: { value: outerRadius },
      outerHeightNegative: { value: outerHeightNegative },
      outerHeightPositive: { value: outerHeightPositive },
      outerYcoord: { value: outerYcoord },
      phaseShift: { value: phaseShift },
      innerRadius: { value: innerRadius },
      innerHeightNegative: { value: innerHeightNegative },
      innerHeightPositive: { value: innerHeightPositive },
      innerYcoord: { value: innerYcoord },
      showInnerCylinder: { value: true },
      showOuterCylinder: { value: true },
      yShift: { value: yShift },
      showCloak: { value: false },
      showLens: { value: true },
      rotAngle: { value: rotAngle },
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
      rectangles: {
        value: rectangles,
      },
      lensSurfaces: { value: lensSurfaces },

      colors: {
        value: [colorRed, colorGreen],
      },

      Camera: {
        value: {
          viewDirection: new THREE.Vector3(0, 0, -1),
          apertureXHat: new THREE.Vector3(1, 0, 0),
          apertureYHat: new THREE.Vector3(0, 1, 0),
          focusDistance: 10.0,
          apertureRadius: infoObject.apertureRadius,
          randomNumbersX: randomNumbersX,
          randomNumbersY: randomNumbersY,
          noOfRays: 1,
        },
      },
      Sphere: {
        value: {
          visible: true,
          centre: new THREE.Vector3(0, 0, 0),
          size: sphereRadius,
        },
      },
    },
    vertexShader: vertexShaderCode,
    fragmentShader: fragmentShaderCode,
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
    rotAngle: rotAngle,
    yShift: yShift,
    outerHeightNegative: outerHeightNegative,
    outerHeightPositive: outerHeightPositive,
    phaseShift: phaseShift,
    innerRadius: innerRadius,
    outerYcoord: outerYcoord,
    innerYcoord: innerYcoord,
    innerHeightNegative: innerHeightNegative,
    innerHeightPositive: innerHeightPositive,
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
    showLens: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.showLens.value =
        !infoObject.raytracingSphereShaderMaterial.uniforms.showLens.value;
      showLensControl.name(
        showLens2String(infoObject.raytracingSphereShaderMaterial)
      );
      console.log("works");
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
    // x1: infoObject.x1,
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
    makeEyeLevel: () => {
      infoObject.raytracingSphereShaderMaterial.uniforms.yShift.value =
        infoObject.camera.position.y;
      console.log(yShift);
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
    .add(GUIParams, "phaseShift", 0, 1, 0.05)
    .name("Hologram Phase shift")
    .onChange((pShift) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.phaseShift.value =
        pShift;
      console.log(pShift);
    });

  resonatorYControl = gui
    .add(GUIParams, "resonatorY", 0, 3, 0.001)
    .name("<i>y</i><sub>cloak</sub>")
    .onChange((y_res) => {
      infoObject.resonatorY = y_res;
      infoObject.raytracingSphereShaderMaterial.uniforms.yShift.value = y_res;
      refreshInfo(infoObject);
      console.log(y_res);
    });

  gui.add(GUIParams, "makeEyeLevel").name("Move resonator to eye level");
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
      lensSurface0.focalLength = r;
    });
  sphereFolder
    .add(GUIParams, "sphereHeight", -1, 1, 0.05)
    .name("<i>h</i><sub>sphere</sub>")
    .onChange((h_sphere) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.sphereHeight.value =
        h_sphere;
      console.log(h_sphere);
    });

  const cloakFolder = gui.addFolder("Axicon Cloak Controls");

  showCloakControl = cloakFolder
    .add(GUIParams, "showCloak")
    .name(showCloak2String(infoObject.raytracingSphereShaderMaterial));

  showOuterCylinderControl = cloakFolder
    .add(GUIParams, "showOuterCylinder")
    .name(showOuterCylinder2String(infoObject.raytracingSphereShaderMaterial));

  cloakFolder
    .add(GUIParams, "outerRadius", 0, 1)
    .name("<i>r</i><sub>outer</sub>")
    .onChange((r_outer) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerRadius.value =
        r_outer;
    });

  cloakFolder
    .add(GUIParams, "outerYcoord", -1, 1)
    .name("<i>y</i><sub>outer</sub>")
    .onChange((y_outer) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.outerYcoord.value =
        y_outer;
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
    .add(GUIParams, "innerYcoord", -1, 1)
    .name("<i>y</i><sub>inner</sub>")
    .onChange((y_inner) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.innerYcoord.value =
        y_inner;
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

  const lensFolder = gui.addFolder("Lens Controls ");

  showLensControl = lensFolder
    .add(GUIParams, "showLens")
    .name(showLens2String(infoObject.raytracingSphereShaderMaterial));

  lensFolder
    .add(GUIParams, "rotAngle", 0, 360, 1)
    .name("Rotation angle")
    .onChange((lens_rot) => {
      infoObject.raytracingSphereShaderMaterial.uniforms.rotAngle.value =
        lens_rot;
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
  gui.add(GUIParams, "No of rays", 1, 100, 1).onChange((n) => {
    infoObject.noOfRays = n;
  });

  if (renderer.xr.enabled) {
    vrControlsVisibleControl = gui.add(GUIParams, "vrControlsVisible");
  }

  // create the GUI mesh at the end to make sure that it includes all controls
  GUIMesh = new HTMLMesh(gui.domElement);
  GUIMesh.visible = false;
  vrControlsVisibleControl.name(guiMeshVisible2String(GUIMesh)); // this can be called only after GUIMesh has been created
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
    document.getElementById("shareButton").src =
      "./textures/icons/shareButtonUnavailable.png";
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
