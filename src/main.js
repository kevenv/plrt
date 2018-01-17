/*
	Precomputed Local Radiance Transfert (PLRT)

	Demo for the computer graphics paper discussion group at McGill University.
	Based on:
			Based on: Precomputed Local Radiance Transfer for Real-Time Lighting Design.
			http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.87.8865&rep=rep1&type=pdf

	by Keven Villeneuve
*/

'use strict'

// Constants
var WINDOW_WIDTH = 768/2;
var WINDOW_HEIGHT = 768/2;

var N_PROBES = 5;
var SEARCH_RADIUS = 0.8;
var LIGHT_POWER = 1.0;
var N_ORDER = 2;
var N_COEFFS = N_ORDER*N_ORDER; //todo: make sure this is right

var N_MONTE_CARLO = 100;
var N_SAMPLES_SH = 50;
var RAY_OFFSET = 1e-18;
var MONTE_CARLO_MAX_BOUNCES = 2;

// Globals
var scene = null;
var camera = null;
var renderer = null;
var controls = null;

var light = {};
var lightChanged = true;
var bounds = {};

var bvh = null;
var triangleObjectMap = {}; // map triangleId to objectId
var trianglesNormal = []; // contains each triangles vertex normals, aligned to bvh output
var objects = [];
var probes = [];
var weights = [];
var PLRTCache = [];
var basisSH = [];

var basicShader = null;
var probeMaterial = new THREE.MeshBasicMaterial( {color: 0x0000ff} );
var probeUnselectedMaterial = new THREE.MeshBasicMaterial( {color: 0x000088} );

var cbox_colors = [new THREE.Vector3(1,1,1),new THREE.Vector3(1,1,1),];

// Events
document.addEventListener("load", onLoad());

function onLoad() {
	//initControls();
	onInit();
}

function onRender() {
	requestAnimationFrame(onRender);
	onUpdate();
	renderer.render(scene, camera);
}

function onInit() {
	// viewport
	renderer = new THREE.WebGLRenderer();
	renderer.setSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	var viewport = document.getElementById("viewport");
	viewport.appendChild(renderer.domElement);

	scene = new THREE.Scene();

	// camera
	camera = new THREE.PerspectiveCamera(37, WINDOW_WIDTH/WINDOW_HEIGHT, 0.1, 1000);
	camera.up.set(0,0,1);
	camera.position.set(0,-4,1);
	camera.lookAt(new THREE.Vector3(0,0,1));

	// controls
	controls = new THREE.OrbitControls(camera, renderer.domElement);
	controls.target.set(0,0,1);

	// shader
	basicShader = createShader();

	onMeshLoad();
	//onAfterLoad();
}

function onMeshLoad() {
	loadCbox();
}

function onAfterLoad() {
	add3DAxis();
	addLight();
	addSearchRadius();
	addProbes();
	//addTestPlane();
	//addTestProbe();
	//addTestOccluder();

	buildBVH(objects);
	//computeBasisSH();
	computeProbes();

	onRender();
}

function addLight() {
	var geometry = new THREE.SphereGeometry(0.2, 32, 32);
	var material = new THREE.MeshBasicMaterial( {color: 0xff0000} );
	light = new THREE.Mesh(geometry, material);
	light.position.z = 0.4;
	//light.position.z = 1.6;
	//light.position.x = -0.3;
	scene.add(light);
}

function addProbe(x,y,z) {
	var geometry = new THREE.SphereGeometry(0.1, 16, 16);
	var material = new THREE.MeshBasicMaterial( {color: 0x0000ff} );
	var probe = new THREE.Mesh(geometry, material); //TODO: use same mesh&material for all probes
	probe.position.set(x,y,z);
	probes.push({
		"mesh": probe, 
		"position": probe.position
	});
	scene.add(probe);
}

function addProbes() {
	var X_W = N_PROBES;
	var Y_W = N_PROBES;
	var Z_W = N_PROBES;
	N_PROBES = X_W*Y_W*Z_W;
	var PROBE_SPACING = 0.5;
	var OFFSET_X = -1;
	var OFFSET_Y = -1;
	var OFFSET_Z = 0.2;
	for(var x = 0; x < X_W; x++) {
		for(var y = 0; y < Y_W; y++) {
			for(var z = 0; z < Z_W; z++) {
				addProbe(OFFSET_X+x*PROBE_SPACING,OFFSET_Y+y*PROBE_SPACING,OFFSET_Z+z*PROBE_SPACING)
			}
		}
	}

	weights = new Array(N_PROBES);
}

function loadCbox() {
	var CBOX_FILE_PATH = "http://localhost:8000/assets/cbox/";
	var cbox_parts = ["cbox_walls","cbox_lwall","cbox_rwall","cbox_bigbox","cbox_smallbox"];
	cbox_colors = [
		new THREE.Color(0.725,0.71,0.68),
		new THREE.Color(0.63,0.065,0.05),
		new THREE.Color(0.14,0.45,0.091),
		new THREE.Color(0.725,0.71,0.68),
		new THREE.Color(0.725,0.71,0.68)
	];

	var loader = new THREE.OBJLoader();

	asyncLoop({
		length : cbox_parts.length,
		functionToLoop : function(loop, i) {
			var filename = CBOX_FILE_PATH + cbox_parts[i] + ".obj";
			loader.load(filename, function(object) {
				var model = object.children[0];

				model.material = basicShader;
				createColorAttrib(model, cbox_colors[i]);

				objects.push(model);
				scene.add(model);
				console.log("loaded " + filename);
				loop();
			});
		},
		callback : function() {
			console.log("add cbox");
			onAfterLoad();
		}
	});
}

function addTestPlane() {
	var geometry = new THREE.PlaneBufferGeometry(2, 2, 10, 10);
	//var material = new THREE.MeshBasicMaterial( {color: 0xffffff, side: THREE.FrontSide} );
	var plane = new THREE.Mesh(geometry, basicShader);
	createColorAttrib(plane, new THREE.Vector3(1.0,1.0,1.0));
	objects.push(plane);
	scene.add(plane);
}

function addTestProbe() {
	var PROBE_SPACING = 0.5;
	var OFFSET_X = -1;
	var OFFSET_Y = -1;
	var OFFSET_Z = 0.2;
	var x=2;
	var y=2;
	var z=2;
	addProbe(OFFSET_X+x*PROBE_SPACING,OFFSET_Y+y*PROBE_SPACING,OFFSET_Z+z*PROBE_SPACING)
}

function addTestOccluder() {
	var geometry = new THREE.BoxBufferGeometry( 0.5, 0.5, 0.05 );
	var cube = new THREE.Mesh( geometry, basicShader);
	createColorAttrib(cube, new THREE.Vector3(1.0,1.0,1.0));

	var rotMat = new THREE.Matrix4();
	//rotMat.makeRotationX(Math.PI/2);

	var verts = cube.geometry.getAttribute("position");
	var N_VERTS = verts.count;
	verts = verts.array;
	for(var v = 0; v < N_VERTS; v++) {
		var vert = new THREE.Vector3(verts[v*3+0], verts[v*3+1], verts[v*3+2]);
		vert.applyMatrix4(rotMat);
		verts[v*3+0] = vert.x;
		verts[v*3+1] = vert.y;
		verts[v*3+2] = vert.z + 1;
	}


	scene.add( cube );
	objects.push(cube);
}

function add3DAxis() {
	var K = 10;
	var matX = new THREE.LineBasicMaterial({color:0xff0000});
	var matY = new THREE.LineBasicMaterial({color:0x00ff00});
	var matZ = new THREE.LineBasicMaterial({color:0x0000ff});
	var geometryX = new THREE.Geometry();
	geometryX.vertices.push(new THREE.Vector3(0, 0, 0));
	geometryX.vertices.push(new THREE.Vector3(K, 0, 0));
	var lineX = new THREE.Line(geometryX, matX);
	scene.add(lineX);
	var geometryY = new THREE.Geometry();
	geometryY.vertices.push(new THREE.Vector3(0, 0, 0));
	geometryY.vertices.push(new THREE.Vector3(0, K, 0));
	var lineY = new THREE.Line(geometryY, matY);
	scene.add(lineY);
	var geometryZ = new THREE.Geometry();
	geometryZ.vertices.push(new THREE.Vector3(0, 0, 0));
	geometryZ.vertices.push(new THREE.Vector3(0, 0, K));
	var lineZ = new THREE.Line(geometryZ, matZ);
	scene.add(lineZ);
}

function addSearchRadius() {
	var boundsG = new THREE.SphereGeometry(SEARCH_RADIUS,8,8);
	var boundsM = new THREE.MeshBasicMaterial( {
		color: new THREE.Color(1,0,1,1),
		wireframe: true
	} );
	bounds = new THREE.Mesh(boundsG, boundsM);
	bounds.position.set(light.position.x,light.position.y,light.position.z);
	scene.add(bounds);
}

function createShader() {
	var vShaderStr = 
	"attribute vec3 mycolor;" +
	"varying vec3 vColor;" + 
	"void main() {" + 
	"	vColor = mycolor;" +
	"	gl_Position = projectionMatrix * modelViewMatrix * vec4(position,1.0);" +
	"}";
	var fShaderStr = 
	"varying vec3 vColor;" +
	"void main() {" + 
	"	gl_FragColor = vec4(vColor,1.0);" +
	"}";
	var shader = new THREE.ShaderMaterial( {
		vertexShader : vShaderStr,
		fragmentShader : fShaderStr
	});
	return shader;
}

function createColorAttrib(mesh, color) {
	var verts = mesh.geometry.getAttribute("position");
	var N_VERTS = verts.count;
	var colors = new Float32Array(N_VERTS * 3);
	for(var i = 0; i < N_VERTS; i++) {
		colors[i*3+0] = color.x;
		colors[i*3+1] = color.y;
		colors[i*3+2] = color.z;
	}
	mesh.geometry.addAttribute("mycolor", new THREE.BufferAttribute(colors, 3));
}



function computeBasisSH() {
	// create SH samples
	basisSH = new Array(N_SAMPLES_SH);
	for(var m = 0; m < N_SAMPLES_SH; m++) {
		basisSH[m] = new Array(N_COEFFS);

		// sample direction
		var sample = new THREE.Vector2(Math.random(), Math.random());
		var w = squareToUniformSphere(sample);
		var pdf = 1.0 / (4.0 * Math.PI);
		var yi = SHEval(w.x, w.y, w.z, N_ORDER);
		for(var i = 0; i < N_COEFFS; i++) {
			basisSH[m][i] = yi[i] / (N_SAMPLES_SH * pdf);
		}
	}
}

function computeProbes() {
	console.log("compute probes...");

	for(var j = 0; j < objects.length; j++) {
		var obj = objects[j];
		var XObj = new Array(N_VERTS); // X matrix / vertex

		var verts = obj.geometry.getAttribute("position");
		var normals = obj.geometry.getAttribute("normal");
		var N_VERTS = verts.count;

		for(var v = 0; v < N_VERTS; v++) {
			
			//alloc X
			var X = new Array(N_COEFFS);
			for(var i = 0; i < N_COEFFS; i++) {
				X[i] = new Array(N_PROBES);
			}

			for(var k = 0; k < N_PROBES; k++) {
				var L_j = MC(v, verts.array, normals.array, probes[k].mesh);

				// project to SH
				//for(var m = 0; m < N_SAMPLES_SH; m++) {
					for(var i = 0; i < N_COEFFS; i++) {
						X[i][k] = L_j;// * basisSH[m][i];
					}
				//}
			}

			var yX = new Array(N_PROBES);
			for(var k = 0; k < N_PROBES; k++) {
				// y^ * X | (1 x nb) * (nb x nl) = dot(y^[i],X[i]) i->N_COEFFS
				var value = 0.0;
				for(var i = 0; i < N_COEFFS; i++) {
					value += X[i][k];//basisSH[0][i] * X[i][k];
				}
				yX[k] = new THREE.Vector3(value,value,value); //todo: color
			}

			XObj[v] = {"X" : X, "yX" : yX, "Lr" : new THREE.Vector3(0,0,0)};
		}

		PLRTCache.push(XObj);
	}

	console.log("[done]");
}

function buildBVH(objects) {
	console.log("build bvh...");

	var triangles = [];

	var triangleIdx = 0;

	for(var i = 0; i < objects.length; i++) {
		var verts = objects[i].geometry.getAttribute("position");
		var normals = objects[i].geometry.getAttribute("normal");
		var N_VERTS = verts.count;
		var verts = verts.array;
		var normals = normals.array;
		for(var k = 0; k < N_VERTS*3; k+=3*3) {
			var v0 = new THREE.Vector3(verts[k+0], verts[k+1], verts[k+2]);
			var v1 = new THREE.Vector3(verts[k+3], verts[k+4], verts[k+5]);
			var v2 = new THREE.Vector3(verts[k+6], verts[k+7], verts[k+8]);
			var n0 = new THREE.Vector3(normals[k+0], normals[k+1], normals[k+2]);
			var n1 = new THREE.Vector3(normals[k+3], normals[k+4], normals[k+5]);
			var n2 = new THREE.Vector3(normals[k+6], normals[k+7], normals[k+8]);
			var triangle = [
				{x: v0.x, y: v0.y, z: v0.z},
				{x: v1.x, y: v1.y, z: v1.z},
				{x: v2.x, y: v2.y, z: v2.z},
			];
			var triangleNormal = [
				{x: n0.x, y: n0.y, z: n0.z},
				{x: n1.x, y: n1.y, z: n1.z},
				{x: n2.x, y: n2.y, z: n2.z},
			];
			triangles.push(triangle);
			trianglesNormal.push(triangleNormal);
			triangleObjectMap[triangleIdx] = i;
			triangleIdx += 1;
		}
	}

	// the maximum number of triangles that can fit in a node before splitting it.
	var maxTrianglesPerNode = 7;
	bvh = new bvhtree.BVH(triangles, maxTrianglesPerNode);

	console.log("[done]");
}

function hitRay(p,w) {
	var its = {};

	var hitObjects = bvh.intersectRay(p, w, true);
	if(hitObjects.length != 0) {
		its.value = true;
		// find out what object was hit
		//todo: sort hit objects
		var objectIdx = triangleObjectMap[hitObjects[0].triangleIndex];
		var obj = objects[objectIdx];
		its.object = obj;

		// hit point
		its.hitPt = hitObjects[0].intersectionPoint;

		// find out normal of hit point
		// todo: vertex normal not face, what's the cloest vertex hit, take this normal?
		var v = 0;
		var n = trianglesNormal[hitObjects[0].triangleIndex][v];
		its.hitNormal = n;
	}

	return its;
}

function squareToUniformSphere(sample) {
	var z = 1.0 - 2.0 * sample.x;
	var r = Math.sqrt(Math.max(0.0, 1.0 - z*z));
	var phi = 2.0 * Math.PI * sample.y;
	return new THREE.Vector3(r * Math.cos(phi), r * Math.sin(phi), z);
}

function MC(v, verts, normals, light) {
	var p = new THREE.Vector3(verts[v*3+0],verts[v*3+1],verts[v*3+2]);
	var n = new THREE.Vector3(normals[v*3+0],normals[v*3+1],normals[v*3+2]);
	// offset ray
	var n_ = n.clone();
	n_.multiplyScalar(RAY_OFFSET);
	p.add(n_);

	var Lr = 0.0;
	for(var i = 0; i < N_MONTE_CARLO; i++) {
		Lr += path_tracer(p,{},n,light,0);
	}
	Lr *= 1.0 / N_MONTE_CARLO;

	return Lr;
}

function path_tracer(woP, woD, n, light, bounces) {
	var hit = hitRay(woP,woD);

	if(bounces > 0 && !hit.value) {
		return 0.0;
	}
	var x = woP;
	if(bounces > 0) {
		n = hit.hitNormal;
		x = hit.hitPt;
	}

	if(bounces >= MONTE_CARLO_MAX_BOUNCES) {
		return 0.0;
	}

	var Lr = 0.0;

	//diffuse (constant added during shading)
	var fr = 1.0;

	//direct
	var wiDir = new THREE.Vector3();
	wiDir.subVectors(light.position, x);
	wiDir.normalize();

	var itsLight = hitRay(x,wiDir);
	var V = !itsLight.value;
	var cosTheta = Math.max(0.0, wiDir.dot(n));
	if(V) {
		var Le = 1.0;
		Lr += fr * Le * cosTheta;
	}

	//indirect
	var sample = new THREE.Vector2(Math.random(), Math.random());
	var wiInd = squareToUniformSphere(sample);
	var pdf = 1.0 / (4.0 * Math.PI);
	//wiInd.add(p);//to world space TODO:
	wiInd.normalize();

	var cosTheta = Math.max(0.0, wiInd.dot(n));

	Lr += (fr * path_tracer(x, wiInd, n, light, ++bounces) * cosTheta) / pdf;

	return Lr;
}

function isProbeVisible(light, probe) {
	var p = light.position;
	var w = new THREE.Vector3();
	w.subVectors(probe.position, light.position);
	var rProbe = w.length();
	w.normalize();

	var its = hitRay(p,w);
	var V = !its.value;
	if(V) return true;

	var rHitVec = new THREE.Vector3();
	rHitVec.subVectors(its.hitPt, light.position);
	var rHit = rHitVec.length();

	return rHit > rProbe; // occluder in front of probe
}

function computeWeights() {
	var sumWeights = 0.0;
	var maxW = 0.0;
	for(var i = 0; i < N_PROBES; i++) {
		probes[i].mesh.material = probeMaterial;

		weights[i] = 0.0;

		var rVec = new THREE.Vector3();
		rVec.subVectors(probes[i].position,light.position)

		var dist = rVec.length();
		if(dist > SEARCH_RADIUS) {
			probes[i].mesh.material = probeUnselectedMaterial;
			continue;
		}

		if(isProbeVisible(light, probes[i])) {
			var w = Math.pow(1.0 - dist/SEARCH_RADIUS, 4);
			weights[i] = w;
			sumWeights += w;
			if(w > maxW) maxW = w;
		} else {
			probes[i].mesh.material = probeUnselectedMaterial;
		}
	}

	for(var i = 0; i < N_PROBES; i++) {
		// normalize
		weights[i] /= sumWeights;
		// multiply by power
		weights[i] *= LIGHT_POWER;
	}

	var sumW=0.0;
	for(var i = 0; i < N_PROBES; i++) {
		sumW += weights[i];
	}
	//console.log(sumW);

	//console.log(weights);
}

function computeVertexRadiance() {
	for(var i = 0; i < objects.length; i++) {
		var obj = objects[i];
		var verts = obj.geometry.getAttribute("mycolor");
		for(var v = 0; v < verts.count; v++) {
			var XObj = PLRTCache[i];
			var yX = XObj[v].yX;
			// yX * w | (1 x nl) * (nl x 1) = dot(yX,w)
			var Lr = new THREE.Vector3(0,0,0);
			var nl = weights.length;
			for(var k = 0; k < nl; k++) {
				var tmp = yX[k].clone();
				tmp.multiplyScalar(weights[k]);
				Lr.add(tmp);
			}

			XObj[v].Lr = Lr;
			//console.log(Lr);
			//XObj[v].Lr = new THREE.Vector3(0,1,0);
		}
	}
}

function updateVertex() {
	for(var i = 0; i < objects.length; i++) {
		var obj = objects[i];
		var verts = obj.geometry.getAttribute("mycolor");
		for(var v = 0; v < verts.count; v++) {

			var XObj = PLRTCache[i];
			var Lr = XObj[v].Lr;

			verts.array[v*3+0] = Lr.x * cbox_colors[i].r/Math.PI;
			verts.array[v*3+1] = Lr.y * cbox_colors[i].g/Math.PI;
			verts.array[v*3+2] = Lr.z * cbox_colors[i].b/Math.PI;
		}

		verts.needsUpdate = true;
	}
}

function moveLight() {
	var time = performance.now() * 0.005;
	light.position.x = Math.sin( time * 0.6 ) * 0.6;
	light.position.y = Math.sin( time * 0.7 ) * 0.6 + 0.1;
	light.position.z = Math.sin( time * 0.8 ) * 0.4 + 0.9;

	bounds.position.x = light.position.x;
	bounds.position.y = light.position.y;
	bounds.position.z = light.position.z;

	lightChanged = true;
}

function onUpdate() {
	controls.update();
	moveLight();
	if(lightChanged) {
		computeWeights();
		computeVertexRadiance();
		updateVertex();	
		lightChanged = false;
	}
}
