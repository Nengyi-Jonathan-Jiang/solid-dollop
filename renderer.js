const renderer = (function(){
    const {gl} = wgllib.fullscreenCanvas(false);
    const {core:{math:{toRad,toDeg,m4},Camera,OrthoCamera,Buffer,VertexArrayObject,Program,Texture,events},createAnimation,gameUtil:{FirstPersonController,CubeMeshGenerator}} = wgllib;
    const {sin,cos,tan,asin,acos,atan,min,max,sqrt,pow,PI,random} = Math;
    events.init();
    // events.init();
    
    gl.clearColor(0,0,0,0);
    
    const shaderProgram = new Program(gl,`
    attribute vec3 a_pos;
    attribute vec3 a_color;
    attribute float a_dark;

    varying vec3 v_color;
    varying float v_dark;
    
    uniform mat4 u_mat;
    
    void main(){
        gl_Position = u_mat * vec4(a_pos, 1.0);
        v_color = a_color;
        v_dark = a_dark;
    }
    `,`
    precision mediump float;
    
    varying vec3 v_color;
    varying float v_dark;

    uniform sampler2D u_texture;
    
    void main() {
        gl_FragColor = vec4((1.0 - v_dark * 1.0) * v_color, 1.0);
    }
    `);
    
    
    var VAO = new VertexArrayObject(gl);
	var VAO2= new VertexArrayObject(gl);
	
    VAO.bind();

    function makeMesh(){
        return new Buffer(gl);
    }
    function bindMesh(VBO){
        VBO.bind();
        VAO.vertexAttribPointer(VBO, shaderProgram.getAttribLoc("a_pos"), "FLOAT", 3, 28, 0);
        VAO.vertexAttribPointer(VBO, shaderProgram.getAttribLoc("a_color"), "FLOAT", 3, 28, 12);
        VAO.vertexAttribPointer(VBO, shaderProgram.getAttribLoc("a_dark"),"FLOAT", 1, 28, 24);
    }

    bindMesh(new Buffer(gl));

	gl.polygonOffset(1,1);
	gl.enable(gl.POLYGON_OFFSET_FILL);
	
    // const camera = new Camera(gl, [0,20,0],[0,0], .01, 1000);
	const camera = new OrthoCamera(gl, [0,20,0],[0,0], .01, -400, 400);
    
    function preDraw(currTime,elapsedTime){
        camera.recompute_projection(toRad(70));
        shaderProgram.uniformMat("u_mat", camera.get_matrix());
    
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    }
	
    function draw(currTime,elapsedTime,numTriangles){
        camera.draw(shaderProgram, VAO, gl.TRIANGLES, 0, numTriangles);
    };

	
    function drawWireFrame(currTime,elapsedTime,numTriangles){
        camera.draw(shaderProgram, VAO, gl.TRIANGLES, 0, numTriangles);
    };

    return {preDraw:preDraw,draw:draw,camera:camera, VAO:VAO,makeMesh:makeMesh,bindMesh:bindMesh};
})();
