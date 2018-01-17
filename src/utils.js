'use strict'

// JSON file read/write
function writeJson(object, name, type) {
    var text = JSON.stringify(object);
    var a = document.createElement("a");
    var file = new Blob([text], {type: type});
    a.href = URL.createObjectURL(file);
    a.download = name;
    a.click();
}

function readJson(file, callback) {
    var rawFile = new XMLHttpRequest();
    rawFile.overrideMimeType("application/json");
    rawFile.open("GET", file, true);
    rawFile.onreadystatechange = function() {
        if (rawFile.readyState === 4) {
            var data = JSON.parse(rawFile.responseText);
            callback(data);
        }
    }
    rawFile.send(null);
}

// Async bullshit
var asyncLoop = function(o)
{
    var i=-1;

    var loop = function() {
        i++;
        if(i==o.length){o.callback(); return;}
        o.functionToLoop(loop, i);
    }
    loop();//init
}

function toRadians(deg) {
    return deg * Math.PI / 180;
};
