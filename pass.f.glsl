#version 300 es
precision highp float;
precision highp int;
in vec3 vcolor;
in vec3 vnormal;
in vec3 wcoord;
in vec3 wvnormal;
in vec3 vdepth;
in vec3 vpos;

in vec3 spot;

uniform float exposure;
uniform vec4 pointlight;
out vec4 color;

void main() {
    vec3 tmp = (normalize(vnormal));
    //uncomment this line for normal map
    //color = vec4(tmp,1);
    //uncomment here for depth map
    //color = vec4(vdepth, 1);
    //point light n*l here
    vec3 l = vec3(pointlight[0] - vpos[0], pointlight[1] - vpos[1], pointlight[2] - vpos[2]);
    vec3 colorrgb = dot(tmp , l)  * vec3(115.0/255.0, 161.0/255.0, 227.0/255.0); //vcolor;
    float r = length(l);
    colorrgb /= r*r*r;
    colorrgb *= pointlight[3];	//pointlight[3] is intensity
    colorrgb = clamp(colorrgb, 0.0, 1.0);
    
    //specular highlight (Phong model)
    float rscale = 2.0 * dot( l / r, tmp);
    vec3 reflect = vec3(0.0, 0.0, 0.0);
    if (rscale > 0.0) {	reflect = rscale * tmp - (l / r);}
    float rdotv = clamp(dot(reflect, - vpos / length(vpos)), 0.0, 1.0);
    float specular = 0.1 * clamp(pow(rdotv , 5.0) * pointlight[3], 0.0, 1.0);
    specular /= r*r;
    colorrgb += specular;
    
    //ambient component
    colorrgb += vec3(0.01, 0.01, 0.01);
    
    colorrgb = clamp(colorrgb, 0.0, 1.0);
    color = vec4(colorrgb, 1);
	
	//uncomment here for mask
	//color = vec4(1, 1, 1, 1);
}
