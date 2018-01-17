'use strict'

function initControls() {
    var text_L_intensity = document.getElementById("text_L_intensity");
    text_L_intensity.value = L_INTENSITY;

    var text_L_direction = document.getElementById("text_L_direction");
    text_L_direction.value = L_ANGLE;
    L_DIR = computeLightDir(L_ANGLE);

    var text_L_r = document.getElementById("text_L_r");
    text_L_r.value = L_r;

    var text_L_d = document.getElementById("text_L_d");
    text_L_d.value = L_d;

    var text_montecarlo = document.getElementById("text_montecarlo");
    text_montecarlo.value = N_MONTE_CARLO;

    var text_order = document.getElementById("text_order");
    text_order.value = N_ORDER;

    var text_savePRT = document.getElementById("text_savePRT");
    text_savePRT.value = PRECOMPUTE_FILE_NAME;

    var sliderIntensity = document.getElementById("slider_L_intensity");
    sliderIntensity.defaultValue = L_INTENSITY;
    sliderIntensity.min = 0.0;
    sliderIntensity.max = 3.0;
    sliderIntensity.step = 0.1;
    sliderIntensity.addEventListener("input", function() {
        L_INTENSITY = parseFloat(sliderIntensity.value);
        var text_L_intensity = document.getElementById("text_L_intensity");
        text_L_intensity.value = L_INTENSITY;
    });

    var sliderDirection = document.getElementById("slider_L_direction");
    sliderDirection.defaultValue = L_ANGLE;
    sliderDirection.min = 0.0;
    sliderDirection.max = 180.0;
    sliderDirection.step = 5.0;
    sliderDirection.addEventListener("input", function() {
        L_ANGLE = parseFloat(sliderDirection.value);
        var text_L_direction = document.getElementById("text_L_direction");
        text_L_direction.value = L_ANGLE;
        L_DIR = computeLightDir(L_ANGLE);
        precomputeL();
    });

    var sliderL_r = document.getElementById("slider_L_r");
    sliderL_r.defaultValue = L_r;
    sliderL_r.min = 0.0;
    sliderL_r.max = 3.0;
    sliderL_r.step = 0.1;
    sliderL_r.addEventListener("input", function() {
        L_r = parseFloat(sliderL_r.value);
        var text_L_r = document.getElementById("text_L_r");
        text_L_r.value = L_r;
        precomputeL();
    });

    var sliderL_d = document.getElementById("slider_L_d");
    sliderL_d.defaultValue = L_d;
    sliderL_d.min = 0.0;
    sliderL_d.max = 3.0;
    sliderL_d.step = 0.1;
    sliderL_d.addEventListener("input", function() {
        L_d = parseFloat(sliderL_d.value);
        var text_L_d = document.getElementById("text_L_d");
        text_L_d.value = L_d;
        precomputeL();
    });

    var button_loadModel = document.getElementById("button_loadModel");
    button_loadModel.addEventListener("click", function() {
        var file_loadModel = document.getElementById("file_loadModel");
        MODEL_FILE_NAME = file_loadModel.value.substring(12,file_loadModel.value.length);
        loadModel();
    });

    var slider_montecarlo = document.getElementById("slider_montecarlo");
    slider_montecarlo.defaultValue = N_MONTE_CARLO;
    slider_montecarlo.min = 0;
    slider_montecarlo.max = 1000;
    slider_montecarlo.step = 100;
    slider_montecarlo.addEventListener("input", function() {
        N_MONTE_CARLO = parseFloat(slider_montecarlo.value);
        var text_montecarlo = document.getElementById("text_montecarlo");
        text_montecarlo.value = N_MONTE_CARLO;
    });

    text_montecarlo.addEventListener("change", function() {
        N_MONTE_CARLO = parseFloat(text_montecarlo.value);
        slider_montecarlo.value = N_MONTE_CARLO;
    });

    var slider_order = document.getElementById("slider_order");
    slider_order.defaultValue = N_ORDER;
    slider_order.min = 3;
    slider_order.max = 5;
    slider_order.step = 1;
    slider_order.addEventListener("input", function() {
        N_ORDER = parseFloat(slider_order.value);
        N_COEFFS = N_ORDER*N_ORDER;
        var text_order = document.getElementById("text_order");
        text_order.value = N_ORDER;
    });

    text_order.addEventListener("change", function() {
        N_ORDER = parseFloat(text_order.value);
        N_COEFFS = N_ORDER*N_ORDER;
        slider_order.value = N_ORDER;
    });

    var button_computePRT = document.getElementById("button_computePRT");
    button_computePRT.addEventListener("click", function() {
        precomputeL();
        precomputeG(false);
    });

    var button_savePRT = document.getElementById("button_savePRT");
    button_savePRT.addEventListener("click", function() {
        var text_savePRT = document.getElementById("text_savePRT");
        PRECOMPUTE_FILE_NAME = text_savePRT.value;
        if(PRTCacheGood) {
            writeJson(PRTCache, PRECOMPUTE_FILE_NAME, 'text/plain');
        }
    });

    var button_loadPRT = document.getElementById("button_loadPRT");
    button_loadPRT.addEventListener("click", function() {
        var file_loadPRT = document.getElementById("file_loadPRT");
        PRECOMPUTE_FILE_NAME = file_loadPRT.value.substring(12,file_loadPRT.value.length);

        var text_savePRT = document.getElementById("text_savePRT");
        text_savePRT.value = PRECOMPUTE_FILE_NAME;
        
        var loc = window.location.pathname;
        var dir = loc.substring(0, loc.lastIndexOf('/'));
        PRECOMPUTE_FILE_PATH = dir + "/" + PRECOMPUTE_FILE_NAME;
        precomputeL();
        precomputeG(true);
    });
}

function computeLightDir(angleDeg) {
    var v = new THREE.Vector3(0,1,0); // at 0 deg
    var rotMat = new THREE.Matrix4();
    rotMat.makeRotationX(toRadians(angleDeg));
    v.applyMatrix4(rotMat);
    return v;
}
