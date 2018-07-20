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
//spotlight input vector

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
    
    //Will move ambient below until after testing/review
    //input S as a vec3 (ASSUME S is normalized)
    //give arbitrary values to a,b,c
    //float a = 1;
    //float b = 1;
    //float c = 1;
    //a, b, c, and distance falloff unnecessary since already calculated for specular and diffuse that are being multiplied
    
    //if angle between L and S (alpha) <= cutoff angle for spotlight (beta)
    //follow through, otherwise, ldotspot == 0
    //remember to normalize l (assume s is normalized)
    //angle between L and S --> cos^-1 (l/r dot spot)/(r* length(spot))
    
    float alpha = acos(dot( l / r, spot)/(r * length(s)));    
    float ldotspot = 0;

    //the power of 5 (currently arbitrary) is the angular falloff coefficient
    if(alpha <= ((3.1415926535897932384626433832795)/2)){
    ldotspot = pow(dot(l, spot), 5);
    }

    //multiply by diffuse and specular
    colorrgb *= ldotspot;

    //ambient component
    colorrgb += vec3(0.01, 0.01, 0.01);
    
    colorrgb = clamp(colorrgb, 0.0, 1.0);
    color = vec4(colorrgb, 1);
	
	//uncomment here for mask
	//color = vec4(1, 1, 1, 1);
}
