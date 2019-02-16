

/*

  #####
 #     # #       ####  #####    ##   #       ####
 #       #      #    # #    #  #  #  #      #
 #  #### #      #    # #####  #    # #       ####
 #     # #      #    # #    # ###### #           #
 #     # #      #    # #    # #    # #      #    #
  #####  ######  ####  #####  #    # ######  ####

*/


var DELTA_TOTAL = 100;
var POINT_JUMP_VALUE = 1;
var DIRECTION_POINT_COUNT = 100;
var ROTATION_JUMP = 1;

var g_globalState = {
    canvasClickLocation: {x: .5, y: .5},
};

var imgw = 200;
var imgh = 200;
var imgsrc = "image8-2.jpg";
var g_img = new Image();
g_img.src = imgsrc;


// #     #
// ##   ##   ##   ##### #    #
// # # # #  #  #    #   #    #
// #  #  # #    #   #   ######
// #     # ######   #   #    #
// #     # #    #   #   #    #
// #     # #    #   #   #    #
//math

//a = [1,0,0], b = [[1],[0],[0]]
//[1,0,0]*[[1],[0],[0]] = [1]
function matrixMultiply(a, b) {
    var aNumRows = a.length, aNumCols = a[0].length,
        bNumRows = b.length, bNumCols = b[0].length,
        m = new Array(aNumRows);  // initialize array of rows
    for (var r = 0; r < aNumRows; ++r) {
        m[r] = new Array(bNumCols); // initialize the current row
        for (var c = 0; c < bNumCols; ++c) {
            m[r][c] = 0;             // initialize the current cell
            for (var i = 0; i < aNumCols; ++i) {
                m[r][c] += a[r][i] * b[i][c];
            }
        }
    }
    return m;
}

function getCosFromDegrees(degrees) {
    return Math.cos(degrees * Math.PI/180);
}

function getSinFromDegrees(degrees) {
    return Math.sin(degrees * Math.PI/180);
}

function inner_bilinear_interp(q00, q10, q01, q11, x, y) {
    var un_x = 1.0 - x; var un_y = 1.0 - y;
    return (q00 * un_x * un_y + q10 * x * un_y + q01 * un_x * y + q11 * x * y);
}

function bilinearInterp(image, xVal, yVal) {
    var x1 = Math.floor(xVal);
    var x2 = Math.floor(xVal) + 1;
    var y1 = Math.floor(yVal);
    var y2 = Math.floor(yVal) + 1;
    if (y2 >= imgh) {
        y2 = imgh - 1;
    }
    if (x2 >= imgw) {
        x2 = imgw - 1;
    }
    return inner_bilinear_interp(image[y1][x1], image[y2][x1], image[y1][x2], image[y2][x2], xVal - x1, yVal - y1)
}


function getScaleMatrix(scaleX, scaleY) {
    return [[scaleX, 0, 0], [0, scaleY, 0], [0, 0, 1]];
}

function getNormScaleMatrix(scale) {
    var normX = Math.sqrt(scale);
    var normY = 1.0 / (Math.sqrt(scale));
    return [[normX, 0, 0], [0, normY, 0], [0, 0, 1]];
}


function getRotationMatrix(inRotation) {
    var toRads = inRotation * Math.PI / 180.0;
    return [
        [Math.cos(toRads), -Math.sin(toRads), 0],
        [Math.sin(toRads), Math.cos(toRads), 0],
        [0, 0, 1]
    ];
}

/*

 #####                                                ######
#     # #####    ##   #####  # ###### #    # #####    #     # ######  ####   ####  ###### #    # #####
#       #    #  #  #  #    # # #      ##   #   #      #     # #      #      #    # #      ##   #   #
#  #### #    # #    # #    # # #####  # #  #   #      #     # #####   ####  #      #####  # #  #   #
#     # #####  ###### #    # # #      #  # #   #      #     # #           # #      #      #  # #   #
#     # #   #  #    # #    # # #      #   ##   #      #     # #      #    # #    # #      #   ##   #
 #####  #    # #    # #####  # ###### #    #   #      ######  ######  ####   ####  ###### #    #   #

*/


function zeros(x) { var r = new Array(x); for (var i = 0; i < x; ++i) { r[i] = 0; } return r; }

function scale(ret, value, c) {
    for (var i = 0; i < value.length; ++i) {
        ret[i] = value[i] * c;
    }
}

function weightedSum(ret, w1, v1, w2, v2) {
    for (var j = 0; j < ret.length; ++j) {
        ret[j] = w1 * v1[j] + w2 * v2[j];
    }
}

/** minimizes a function using the downhill simplex method */
function nelderMead(f, x0, parameters) {
    parameters = parameters || {};

    var maxIterations = parameters.maxIterations || x0.length * 200,
        nonZeroDelta = parameters.nonZeroDelta || 1.05,
        zeroDelta = parameters.zeroDelta || 0.001,
        minErrorDelta = parameters.minErrorDelta || 1e-6,
        minTolerance = parameters.minErrorDelta || 1e-5,
        rho = (parameters.rho !== undefined) ? parameters.rho : 1,
        chi = (parameters.chi !== undefined) ? parameters.chi : 2,
        psi = (parameters.psi !== undefined) ? parameters.psi : -0.5,
        sigma = (parameters.sigma !== undefined) ? parameters.sigma : 0.5,
        maxDiff;

    // initialize simplex.
    var N = x0.length,
        simplex = new Array(N + 1);
    simplex[0] = x0;
    simplex[0].fx = f(x0);
    simplex[0].id = 0;
    for (var i = 0; i < N; ++i) {
        var point = x0.slice();
        point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
        simplex[i+1] = point;
        simplex[i+1].fx = f(point);
        simplex[i+1].id = i+1;
    }

    function updateSimplex(value) {
        for (var i = 0; i < value.length; i++) {
            simplex[N][i] = value[i];
        }
        simplex[N].fx = value.fx;
    }

    var sortOrder = function(a, b) { return a.fx - b.fx; };

    var centroid = x0.slice(),
        reflected = x0.slice(),
        contracted = x0.slice(),
        expanded = x0.slice();

    for (var iteration = 0; iteration < maxIterations; ++iteration) {
        simplex.sort(sortOrder);

        if (parameters.history) {
            // copy the simplex (since later iterations will mutate) and
            // sort it to have a consistent order between iterations
            var sortedSimplex = simplex.map(function (x) {
                var state = x.slice();
                state.fx = x.fx;
                state.id = x.id;
                return state;
            });
            sortedSimplex.sort(function(a,b) { return a.id - b.id; });

            parameters.history.push({x: simplex[0].slice(),
                fx: simplex[0].fx,
                simplex: sortedSimplex});
        }

        maxDiff = 0;
        for (i = 0; i < N; ++i) {
            maxDiff = Math.max(maxDiff, Math.abs(simplex[0][i] - simplex[1][i]));
        }

        if ((Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) &&
            (maxDiff < minTolerance)) {
            break;
        }

        // compute the centroid of all but the worst point in the simplex
        for (i = 0; i < N; ++i) {
            centroid[i] = 0;
            for (var j = 0; j < N; ++j) {
                centroid[i] += simplex[j][i];
            }
            centroid[i] /= N;
        }

        // reflect the worst point past the centroid  and compute loss at reflected
        // point
        var worst = simplex[N];
        weightedSum(reflected, 1+rho, centroid, -rho, worst);
        reflected.fx = f(reflected);

        // if the reflected point is the best seen, then possibly expand
        if (reflected.fx < simplex[0].fx) {
            weightedSum(expanded, 1+chi, centroid, -chi, worst);
            expanded.fx = f(expanded);
            if (expanded.fx < reflected.fx) {
                updateSimplex(expanded);
            }  else {
                updateSimplex(reflected);
            }
        }

        // if the reflected point is worse than the second worst, we need to
        // contract
        else if (reflected.fx >= simplex[N-1].fx) {
            var shouldReduce = false;

            if (reflected.fx > worst.fx) {
                // do an inside contraction
                weightedSum(contracted, 1+psi, centroid, -psi, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < worst.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            } else {
                // do an outside contraction
                weightedSum(contracted, 1-psi * rho, centroid, psi*rho, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < reflected.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            }

            if (shouldReduce) {
                // if we don't contract here, we're done
                if (sigma >= 1) break;

                // do a reduction
                for (i = 1; i < simplex.length; ++i) {
                    weightedSum(simplex[i], 1 - sigma, simplex[0], sigma, simplex[i]);
                    simplex[i].fx = f(simplex[i]);
                }
            }
        } else {
            updateSimplex(reflected);
        }
    }

    simplex.sort(sortOrder);
    return {fx : simplex[0].fx,
        x : simplex[0]};
}

function pythagoras(xval, yval) {
    return Math.sqrt((xval**2)+(yval**2));
}

function getSquaredDistanceOfPoint(xval, yval) {
     return pythagoras(xval, yval) ** 2;
}

function getTotalDiffSquaredOfPoints_matrix(points) {
    var total = 0;
    for (var i = 0; i < points.length; i++) {
        total += getSquaredDistanceOfPoint(points[i][0], points[i][1]);
    }
    return total;
}

function rotatePoint_matrix(degrees, point) {
    rads = degrees  * Math.PI / 180.0; //convert to rads
    sinT = Math.sin(rads);
    cosT = Math.cos(rads);

    rotMat = [[cosT,sinT],[-sinT,cosT]];
    pointMat = [[point[0]], [point[1]]];

    return matrixMultiply(rotMat, pointMat);
}

function applyTransformToPoint_matrix(degrees, normX, point) {
    var ret = point;
    ret = rotatePoint_matrix(degrees, ret);

    ret = [ ret[0]*normX, ret[1] ];

    ret = rotatePoint_matrix(-degrees, ret);
    return ret
}

function applyTransformToAllPoints(tetha, normX, normY, points) {
    var ret = [];

    for (var i = 0; i < points.length; i++){

        newPoint = points[i];
        newPoint = applyTransformToPoint_matrix(tetha, normX, [newPoint[0], newPoint[1]]);
        newPoint = applyTransformToPoint_matrix(tetha+90, normY, [newPoint[0], newPoint[1]]);

        ret.push( [ newPoint[0][0], newPoint[1][0] ] )
    }
    return ret
}

function getCenterPointOfPoly(arr) {
    var minX, maxX, minY, maxY;
    for(var i=0; i< arr.length; i++){
        minX = (arr[i][0] < minX || minX == null) ? arr[i][0] : minX;
        maxX = (arr[i][0] > maxX || maxX == null) ? arr[i][0] : maxX;
        minY = (arr[i][1] < minY || minY == null) ? arr[i][1] : minY;
        maxY = (arr[i][1] > maxY || maxY == null) ? arr[i][1] : maxY;
    }
    return [(minX + maxX) /2, (minY + maxY) /2];
}

function getTranslateMatrix(x, y) {
    return [
        [1, 0, x],
        [0, 1, y],
        [0, 0, 1]
    ];
}

function pointsToMatrixPoints(shape) {
    var outShape = [];
    for (var i = 0; i < shape[0].length; i++){
        var line = [];
        for (var j = 0; j < shape.length; j++) {
            line.push(shape[j][i])
        }
        outShape.push(line);
    }
    var line = [];
    for (var j = 0; j < shape.length; j++) {
        line.push(1)
    }
    outShape.push(line);

    return outShape;
}

function matrixPointsToPoints(shape) {
    var outShape = [];
    for (var i = 0; i < shape[0].length; i++){
        outShape.push( [ shape[0][i], shape[1][i] ] );
    }
    return outShape;
}

function prepForMinimise(shape) {
    var trans = getCenterPointOfPoly(shape);
    var transmat = getTranslateMatrix(-trans[0], -trans[1]);
    var points_1 = pointsToMatrixPoints(shape);
    var points_2 = matrixMultiply(transmat, points_1);
    return matrixPointsToPoints(points_2);
}

function scaleInDirection(shape, angle, scale) {
    var normX = Math.sqrt(scale);
    var normY = 1.0 / (Math.sqrt(scale));
    return applyTransformToAllPoints(angle, normX, normY, shape)
}

function getAngleForOnePoint_matrix(point) {

    if(point[0] === 0 && point[1] >= 0) {
        return 270;
    } else if(point[0] === 0 && point[1] < 0) {
        return 90;
    }

    const atanVal = Math.atan(point[1]/point[0]);
    let degs = Math.abs(atanVal * 180.0/Math.PI);

    if (point[1] >= 0 && point[0] >= 0) {
        degs = 360 - degs;
    } else if (point[1] < 0 && point[0] >= 0) {
        //degs = degs;
    } else if (point[1] >= 0 && point[0] < 0) {
        degs += 180;
    } else if (point[1] < 0 && point[0] < 0) {
        degs = 180 - degs;
    }

    return degs
}

function getAngleBetweenTwoPoints_matrix(point1, point2) {
    return Math.abs(getAngleForOnePoint_matrix(point1) - getAngleForOnePoint_matrix(point2))
}

function getTotalDiffSquaredOfAngles_matrix(shape) {
    var retTotal = 0;
    for (var i = 0; i < shape.length; i++) {
        var allAnglesForPoint = [];
        for (var j = 0; j < shape.length; j++) {
            if (j === i)
                continue;

            var angleBetweenPoints = getAngleBetweenTwoPoints_matrix(shape[i], shape[j])
            allAnglesForPoint.push(angleBetweenPoints);
        }
        var lowest = allAnglesForPoint[0];
        for (var k = 1; k < allAnglesForPoint.length; k++) {
            if (allAnglesForPoint[k] < lowest) {
                lowest = allAnglesForPoint[k];
            }
        }
        retTotal += lowest*lowest//should we square this????
    }
    return retTotal
}

function calcDiffSquaredOfEveryPoint(shape, angle, scale) {
    if(scale < 0){
        scale = scale*-1
    }

    newShape = scaleInDirection(shape, angle, scale);
    totalDiff = getTotalDiffSquaredOfPoints_matrix(newShape);
    totalDiff += getTotalDiffSquaredOfAngles_matrix(newShape);
    return totalDiff
}



var g_global_shape = prepForMinimise([[0,0],[0,2],[2,2],[2,0]]);

function scaleAndAngleFunc(X) {
    var shape = g_global_shape;
    var angle = X[0];
    var scale = X[1];
    return calcDiffSquaredOfEveryPoint(shape, angle, scale);
};

//public api
function getValuesToNormaliseScale1(shape) {
    g_global_shape = prepForMinimise(shape);
    let solution = nelderMead(scaleAndAngleFunc, [10, 1]);

    console.log("solution is at " + solution.x);

    return {angle: solution.x[0], scale: Math.abs(solution.x[1])};
}


// #     #                         ###
// #     #  ####  ###### #####      #  #    # #####  #    # #####
// #     # #      #      #    #     #  ##   # #    # #    #   #
// #     #  ####  #####  #    #     #  # #  # #    # #    #   #
// #     #      # #      #####      #  #  # # #####  #    #   #
// #     # #    # #      #   #      #  #   ## #      #    #   #
//  #####   ####  ###### #    #    ### #    # #       ####    #
//user input



function getCurrentPageMousePosition(e) {
    return {
        x: e.pageX,
        y: e.pageY
    };
}


function getCurrentCanvasMousePosition(e, canvasElem) {
    if (e.originalEvent.changedTouches != null && canvasElem != null) {
        var rect = canvasElem.getBoundingClientRect();
        return {
            x: e.originalEvent.changedTouches[0].clientX - rect.left,
            y: e.originalEvent.changedTouches[0].clientY - rect.top
        };
    } else if (e.clientX || e.clientX === 0 && canvasElem != null) {
        var rect = canvasElem.getBoundingClientRect();
        return {
            x: e.clientX - rect.left,
            y: e.clientY - rect.top
        };
    } else {
        console.log("Error: Invalid state");
    }
}

$("#" + "canvasImg1UI").mousedown(function (e) {
    e.preventDefault();

    var canvasElem = $("#" + "canvasImg1")[0];
    const pageMousePosition = getCurrentPageMousePosition(e);
    const canvasMousePosition = getCurrentCanvasMousePosition(e, canvasElem);
    console.log(pageMousePosition);
    console.log(canvasMousePosition);

    g_globalState.canvasClickLocation = {x: canvasMousePosition.x/imgw, y: canvasMousePosition.y/imgh};
    draw()
});

$("#" + "canvasImg2UI").mousedown(function (e) {
    e.preventDefault();

    var canvasElem = $("#" + "canvasImg2")[0];
    const pageMousePosition = getCurrentPageMousePosition(e);
    const canvasMousePosition = getCurrentCanvasMousePosition(e, canvasElem);
    console.log(pageMousePosition);
    console.log(canvasMousePosition);

    g_globalState.canvasClickLocation = {x: canvasMousePosition.x/imgw, y: canvasMousePosition.y/imgh};
    draw()
});


function tomdrawlin() {
    var c = document.getElementById("myCanvas");
    var ctx = c.getContext("2d");
    ctx.beginPath();
    ctx.moveTo(0, 0);
    ctx.lineTo(300, 150);
    ctx.stroke();
}


// #####
// #     # #####    ##   #    #
// #     # #    #  #  #  #    #
// #     # #    # #    # #    #
// #     # #####  ###### # ## #
// #     # #   #  #    # ##  ##
// #####   #    # #    # #    #
//draw


function getPixels(image, width, height, rot, scale, scaleRot) {
    var canvas = document.getElementById('imageMods');
    var ctx = canvas.getContext('2d');
    ctx.save();

    ctx.translate(imgw/2, imgh/2);
    ctx.rotate(scaleRot * Math.PI / 180);

    var normX = Math.sqrt(scale);
    var normY = 1.0 / (Math.sqrt(scale));
    ctx.scale(normX, normY);

    ctx.rotate(-scaleRot * Math.PI / 180);
    ctx.translate(-imgw/2, -imgh/2);

    ctx.translate(imgw/2, imgh/2);
    ctx.rotate(rot * Math.PI / 180);
    ctx.translate(-imgw/2, -imgh/2);
    ctx.drawImage(image, 0, 0, imgw, imgh);
    ctx.restore();
    return ctx.getImageData(0, 0, width, height).data;
}

function toBlackAndWhite(imageData) {
    var output = []

    for (i = 0; i < imgh; i++) {
        var arr = []
        for (j = 0; j < imgw; j++) {
            var index = ( i * (imgw * 4) ) + (j * 4)
            var val = ((imageData[index] + imageData[index+1] + imageData[index+2])/3.0)
            arr.push(val)
        }
        output.push(arr)
    }
    return output
}

function main() {
    //g_p = toBlackAndWhite(getPixels())
    draw();
}

function plotPoints(inputArr, yval, xval) {
    labels = [];
    zvals = [];
    nulllist = [];
    nulllistr = [];
    nulllistl = [];
    var i,j;

    for (i = 0; i < inputArr[0].length; i++)
    {
        labels.push(i);
        zvals.push(inputArr[0][i].z)

        //nulllist.push( (i != xval)? null : inputArr[xval] );
        //nulllistr.push( (i != rightVal)? null : inputArr[rightVal] );
        //nulllistl.push( (i != leftVal)? null : inputArr[leftVal] );
    }
    for (j = 0; j < inputArr[1].length; j++) {
        labels.push(i+ j);
        zvals.push(inputArr[1][j].z)
    }
    // Initialize a Line chart in the container with the ID chart1
    new Chartist.Line('#chart1', {
        labels: labels,
        series: [
            zvals
            //    nulllist,
            //    nulllistr,
            //    nulllistl
        ]
    });
}

function findRightPointOfAnchor(output) {
    var deltaTotal = 0;
    for (i = 0; i < output.length - 1; i++) {
        deltaTotal += Math.abs( output[i+1].z - output[i].z );
        if (deltaTotal > DELTA_TOTAL){
            return output[i];
        }
    }
}

function getDirectionPointsWithjump(xval, yval, xjump, yjump) {
    var dx = xval;
    var dy = yval;
    var output = [];
    var count = 0;
    while(dx < imgw && dx > 0 && dy < imgh && dy > 0) {
        if (count > DIRECTION_POINT_COUNT) {
            break;
        }
        output.push({x: dx, y: dy});
        count++;
        dx = dx + xjump;
        dy = dy + yjump;
    }
    return output;
}

function to_matrix_shape(shape) {
    var ret = [];
    for (let i = 0; i < shape.length; i++) {
        ret.push([shape[i].x, shape[i].y]);
    }
    return ret;
}

function getDirectionPoints(xval, yval, rot) {
    var jumpH = POINT_JUMP_VALUE;
    var cosval = getCosFromDegrees(rot);
    var sinval = getSinFromDegrees(rot);
    return getDirectionPointsWithjump(
        xval + cosval, yval + sinval, cosval*jumpH, sinval*jumpH)
}

function getZValues(image, xval, yval, rot) {
    var points = getDirectionPoints(xval, yval, rot);
    var ret = [];
    for (var i = 0; i < points.length; i++) {
        ret.push( {
            x: points[i].x,
            y: points[i].y,
            z: bilinearInterp(image, points[i].x, points[i].y)
        } );
    }
    return ret;
}

function drawPoint(interactiveCanvasContext, point, colour) {
    interactiveCanvasContext.beginPath();
    interactiveCanvasContext.strokeStyle = colour;
    interactiveCanvasContext.rect(point.x, point.y, 3, 3);
    interactiveCanvasContext.closePath();
    interactiveCanvasContext.stroke();
}

function plotImagePoints(output, xval, yval) {

    plotPoints( output, yval, xval )
}

function drawImageWithTransformations(ctx, img, imgw, imgh, scale, scaleRot) {

    ctx.save();

    ctx.translate(imgw/2, imgh/2);
    ctx.rotate(scaleRot * Math.PI / 180);

    var normX = Math.sqrt(scale);
    var normY = 1.0 / (Math.sqrt(scale));
    ctx.scale(normX, normY);

    ctx.rotate(-scaleRot * Math.PI / 180);
    ctx.translate(-imgw/2, -imgh/2);

    //ctx.translate(imgw/2, imgh/2);
    //ctx.rotate(45 * Math.PI / 180);
    //ctx.translate(-imgw/2, -imgh/2);
    //ctx.translate(200, -200);
    ctx.drawImage(img, 0, 0, imgw, imgh);
    ctx.restore();
}

function getHitPoints(img, imageData, m_xval, m_yval) {

    var blackandwhite = toBlackAndWhite(imageData);
    var rot = 0;
    var shape = [];
    for (let rot = 0; rot < 360; rot += ROTATION_JUMP) {

        var zvals = getZValues(blackandwhite, m_xval, m_yval, rot);
        var hitPoint = findRightPointOfAnchor(zvals);
        if (hitPoint != undefined) {
            shape.push(hitPoint);
        }

    }
    return shape;
}

function getTransformationMatrixFromScale(scale, rotation, imgw, imgh) {
    var mat = getTranslateMatrix(imgw/2, imgh/2);
    mat = matrixMultiply(mat, getRotationMatrix(rotation));
    mat = matrixMultiply(mat, getNormScaleMatrix(scale));
    mat = matrixMultiply(mat, getRotationMatrix(-rotation));
    mat = matrixMultiply(mat, getTranslateMatrix(-imgw/2, -imgh/2));
    return mat;
}

function applyTransformationMatrixToPoint(point, mat) {
    var resPoint = matrixMultiply( mat, [[point.x], [point.y], [1]]);
    return { x: resPoint[0][0], y: resPoint[1][0]};
}


function getTransformedPoint(point, undoScale, undoScaleRotation, scale, scaleRotation, imgw, imgh) {
    var mat = getTransformationMatrixFromScale( undoScale, -undoScaleRotation, imgw, imgh );
    var mat2 = getTransformationMatrixFromScale( scale, scaleRotation, imgw, imgh );
    mat = matrixMultiply(mat, mat2);
    return applyTransformationMatrixToPoint(point, mat);
}

function draw() {

    var scaleRotation1 = -45;
    var scale1 = 2;

    var scaleRotation2 = 45;
    var scale2 = 2;

    var c = document.getElementById("canvasImg1");
    var ctx = c.getContext("2d");
    ctx.clearRect(0, 0, c.width, c.height);
    drawImageWithTransformations(ctx, g_img, imgw, imgh, scale1, scaleRotation1);

    var c2 = document.getElementById("canvasImg2");
    var ctx2 = c2.getContext("2d");
    ctx2.clearRect(0, 0, c2.width, c2.height);
    drawImageWithTransformations(ctx2, g_img, imgw, imgh, scale2, scaleRotation2);

    var cUI = document.getElementById("canvasImg1UI");
    var ctxUI = cUI.getContext("2d");
    ctxUI.clearRect(0, 0, c2.width, c2.height);

    var c2UI = document.getElementById("canvasImg2UI");
    var ctx2UI = c2UI.getContext("2d");
    ctx2UI.clearRect(0, 0, c2.width, c2.height);

    var m_xval = g_globalState.canvasClickLocation.x*imgw;
    var m_yval = g_globalState.canvasClickLocation.y*imgh;
    var m_yval = Math.round(m_yval);
    var m_xval = Math.round(m_xval);

    //var imageData = getPixels(img, imgw, imgh, 120, 2.5);
    var img1HitPoints = getHitPoints(c, ctx.getImageData(0, 0, c.width, c.height).data,
        m_xval, m_yval);

    drawPoint(ctxUI, {x: m_xval, y: m_yval}, "blue");
    for (let i = 0; i < img1HitPoints.length; i++) {
        drawPoint(ctxUI, img1HitPoints[i], "red")
    }

    var transformedPoint = getTransformedPoint({x: m_xval, y: m_yval}, scale1,
                    scaleRotation1, scale2, scaleRotation2, imgw, imgh);

    //var imageData = getPixels(img, imgw, imgh, 120, 2.5);
    var img2HitPoints = getHitPoints(c2, ctx2.getImageData(0, 0, c2.width, c2.height).data,
        transformedPoint.x, transformedPoint.y);

    //FIXME: we can't draw points on the canvas we test
    drawPoint(ctx2UI, transformedPoint, "blue");

    for (let i = 0; i < img2HitPoints.length; i++) {
        drawPoint(ctx2UI, img2HitPoints[i], "red")
    }

    var c1Result = document.getElementById("canvasImg1Result");
    var ctx1Result = c1Result.getContext("2d");
    ctx1Result.clearRect(0, 0, c2.width, c2.height);

    var c2Result = document.getElementById("canvasImg2Result");
    var ctx2Result = c2Result.getContext("2d");
    ctx2Result.clearRect(0, 0, c2.width, c2.height);

    var canvasOverlay = document.getElementById("canvasOverlay");
    var ctxCanvasOverlay = canvasOverlay.getContext("2d");
    ctxCanvasOverlay.clearRect(0, 0, c2.width, c2.height);

    if (img1HitPoints.length > 3) {
        var normVals = getValuesToNormaliseScale1(to_matrix_shape(img1HitPoints));
        drawImageWithTransformations(ctx1Result, c, imgw, imgh, normVals.scale, normVals.angle);
        ctxCanvasOverlay.globalAlpha = 0.5;
        drawImageWithTransformations(ctxCanvasOverlay, c, imgw, imgh, normVals.scale, normVals.angle);
    }

    if (img2HitPoints.length > 3) {
        var normVals = getValuesToNormaliseScale1(to_matrix_shape(img2HitPoints));
        drawImageWithTransformations(ctx2Result, c2, imgw, imgh, normVals.scale, normVals.angle);
        ctxCanvasOverlay.globalAlpha = 0.5;
        drawImageWithTransformations(ctxCanvasOverlay, c2, imgw, imgh, normVals.scale, normVals.angle);
    }


}



