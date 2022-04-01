class Block{
    /** @param {number} id @param {number} light @param {number} data */
    constructor(id,light,data=0){
        this.id = id;
        this.data = data;
        this.light = light;
    }

    static DEBUG = 0;
    static GRASS = 1;
    static DIRT = 2;
    static STONE = 3;

    static getTransparency(blockId){
        return ({    //This object stores how much light passes through the block
            0: 1,
            28: 1,
        }[blockId] || 0);
    }

    static DIFFUSION = 32;
}

class Chunk{
    // static WIDTH = 64; 
    // static HEIGHT = 256;
    // static DEPTH = 64;
	static WIDTH = 8; 
    static HEIGHT = 16;
    static DEPTH = 8;
    static TOTAL_BLOCKS = Chunk.WIDTH * Chunk.HEIGHT * Chunk.DEPTH;

    static MESH_GEN = new wgllib.gameUtil.CubeMeshGenerator(16, 16);

    /** @param {[number,number,number,number][]} [positions] */
    constructor(positions){
		this.blockExists =  new Uint8Array(Chunk.TOTAL_BLOCKS);
		
		this.blockColorR =  new Uint8Array(Chunk.TOTAL_BLOCKS);
		this.blockColorG =  new Uint8Array(Chunk.TOTAL_BLOCKS);
		this.blockColorB =  new Uint8Array(Chunk.TOTAL_BLOCKS);

        this.setFromPositions(positions||[]);
    }

    /** @param {[number,number,number,number][]} positions */
    setFromPositions(positions){
		this.blockExists.fill(0);
		this.blockColorR.fill(0);
		this.blockColorG.fill(0);
		this.blockColorB.fill(0);

        for(let [x,y,z,r,g,b] of positions){
			this.setBlock(x,y,z,[r,g,b]);
		}
    }

	setBlock(x,y,z,[r,g,b]){
		if(!this._inRange(x,y,z)) return;
		this.setBlockExistsAt(x,y,z,1);
		this.setBlockColorAt(x,y,z,[r,g,b]);
	}
	breakBlock(x,y,z){
		if(!this._inRange(x,y,z)) return;
		this.setBlockExistsAt(x,y,z,0);
	}
	
    getMesh(tx,ty,tz){
        const {m4,m4:{dot,add,cross},toRad} = wgllib.core.math;
        const {sin,cos,round} = Math;

        let {WIDTH,HEIGHT,DEPTH,MESH_GEN} = Chunk;

        //Count faces to be rendered (so as not to allocate extra memory)
        let faceCount = 0;
        //For each block
        for(let x = 0; x < WIDTH; x++){
            for(let y = 0; y < HEIGHT; y++){
                for(let z = 0; z < DEPTH; z++){
                    //If it's air (nothing is drawn), don't do anything
                    if(this.getBlockExistsAt(x,y,z) == 0) continue;
                    //for each face
                    for(let [dx,dy,dz] of [[0,1,0],[0,-1,0],[0,0,-1],[0,0,1],[1,0,0],[-1,0,0]]){
                        //get position of adjacent block
                        let [xx,yy,zz] = [x + dx, y + dy, z + dz];
                        if(!this._inRange(xx,yy,zz) || !this.getBlockExistsAt(xx,yy,zz)){
							//Increase faceCount
							faceCount++;
                        }
                    }
                }
            }
        }

        //Normals to faces
        const fNormals = [[0,1,0],[0,-1,0],[0,0,-1],[0,0,1],[1,0,0],[-1,0,0]];
        
        //Initialize array
        let data = new Float32Array(faceCount * 42);
        //Index at which to insert data
        let i = 0;
        //For each block
        for(let x = 0; x < WIDTH; x++){
            for(let y = 0; y < HEIGHT; y++){
                for(let z = 0; z < DEPTH; z++){
                    //If the block is air, do nothing
                    if(this.getBlockExistsAt(x,y,z) == 0) continue;
					
                    //Block color
                    const blockColor = this.getBlockColorAt(x,y,z);
                    
                    //For each face
                    for(let face = 0; face < 6; face++){
                        //Normal vector of that face
                        let faceNormal = fNormals[face];
                        let [dx,dy,dz] = faceNormal;
                        //Block adjacent to that face
                        let [xx,yy,zz] = [x + dx, y + dy, z + dz];
                        if(!this._inRange(xx,yy,zz) || !this.getBlockExistsAt(xx,yy,zz)){
                            //For each vertex of the face
                            for(let vertex = 0; vertex < 6; vertex++){
                                //Set the data (position)
                                [data[i++],data[i++],data[i++]] = MESH_GEN.getPos(face,vertex,[x - tx,y - ty,z - tz]);
                                //Set the data (color)
                                [data[i++],data[i++], data[i++]] = blockColor.map(i=>i/255);
                                //Set the data (light level)
                                data[i++] = 1 - [1,.64,.7,.7,.8,.8][face];
                            }
                        }
                    }
                }
            }
        }

        return data;
    }

    logLighting(){
        function f(n){return " .:-=+*Q#M%@░▒▓█"[Math.floor(n / 16)].repeat(2)}
        for(let y = Chunk.HEIGHT; y --> 0;)
            console.log(
                "%c" + new Array(Chunk.WIDTH).fill(0).map(
                    (_,x)=>new Array(Chunk.DEPTH).fill(0).map(
                        (_,z)=>this.getBlockLightAt(x,y,z)
                    ).map(f).join("")
                ).join("\n"),
                "font-family:monospace;"
            )
    }

    _mapPos(x,y,z){return (x * Chunk.DEPTH + z) * Chunk.HEIGHT + y}
    _inRange(x,y,z){return x>=0&&y>=0&&z>=0&&x<Chunk.WIDTH&&y<Chunk.HEIGHT&&z<Chunk.DEPTH}

	getBlockExistsAt(x,y,z){return this.blockExists[this._mapPos(x,y,z)]}
    setBlockExistsAt(x, y, z, exists) {this.blockExists[this._mapPos(x,y,z)] = exists}
	
	getBlockColorRAt(x,y,z){return this.blockColorR[this._mapPos(x,y,z)]}
    setBlockColorRAt(x,y,z,v) {this.blockColorR[this._mapPos(x,y,z)] = v}
	getBlockColorGAt(x,y,z){return this.blockColorG[this._mapPos(x,y,z)]}
    setBlockColorGAt(x,y,z,v) {this.blockColorG[this._mapPos(x,y,z)] = v}
	getBlockColorBAt(x,y,z){return this.blockColorB[this._mapPos(x,y,z)]}
    setBlockColorBAt(x,y,z,v) {this.blockColorB[this._mapPos(x,y,z)] = v}
	getBlockColorAt(x,y,z){return [this.getBlockColorRAt(x,y,z),this.getBlockColorGAt(x,y,z),this.getBlockColorBAt(x,y,z)]}
	setBlockColorAt(x,y,z,[r,g,b]){this.setBlockColorRAt(x,y,z,r),this.setBlockColorGAt(x,y,z,g),this.setBlockColorBAt(x,y,z,b)}
}