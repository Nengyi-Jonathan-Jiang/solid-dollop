const {core:{math:{toDeg,toRad}}} = wgllib;

var camera = renderer.camera;
[camera.position,camera.rotation] = [[Chunk.WIDTH / 2, 1, Chunk.DEPTH / 2], [0.6154797086703873, 0.7853981633974483]];
camera.scale = .01;

window.onwheel = e=>{
	let dy = e.deltaY / 100;
	camera.scale *= (1 + dy * .1);
}

var VBO = renderer.makeMesh();

let seed = 264;
noise.seed(seed);
let positions = [];
{
    let scale = 12, offset = 4, yScale = 4;
    for(let x = 0; x < Chunk.WIDTH; x++){
        for(let z = 0; z < Chunk.DEPTH; z++){
            // let height = ~~(noise.simplex2(x / scale, z / scale) * yScale / 2 + yScale / 2 + offset);
            // for(let y = 0; y < height; y++) positions.push([x,y,z,255,0,0]);
            // positions.push([x,height,z,255,195,63]);
			if((x + z) % 2 == 0)
				positions.push([x,0,z,225,225,225]);
			else
				positions.push([x,0,z,160,160,160]);
        }
    }
}

let c = new Chunk(positions);
let c_mesh = c.getMesh(0,0,0);

var control = new wgllib.gameUtil.FirstPersonController(renderer.camera,0);
control.flySpeed = 8;

wgllib.createAnimation(function(currTime,elapsedTime){
    control.update(elapsedTime);
    renderer.preDraw();
    renderer.bindMesh(VBO);
	VBO.setData(c_mesh);
	renderer.draw(currTime,elapsedTime, VBO.bytes / 28);
});

window.onmousedown = e=>{
	const {sin,cos} = Math;
	let rot = camera.rot;
	let [a,b] = rot;
	let [x,y,z] = [0,0,1];
	[y,z] = [y * cos(a) + z * sin(a), z * cos(a) - y * sin(a)];
	[x,z] = [x * cos(b) - z * sin(b), z * cos(b) + x * sin(b)];
	console.log(x,y,z);

	let placePos = [x * 3,y * 3,z * 3];
	let truePos = [~~(placePos[0] + Chunk.WIDTH / 2), ~~(placePos[1] + Chunk.HEIGHT / 2), ~~(placePos[2] + Chunk.DEPTH / 2)];

	console.log(truePos);

	switch(e.button){
		case 2:
			c.setBlock(...truePos, [255,255,255]);
			break;
		case 0:
			c.breakBlock(...truePos);
	}
	
	c_mesh = c.getMesh(0,0,0);
}

window.oncontextmenu = e=>e.preventDefault();