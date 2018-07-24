precision highp float;
precision highp int;

//vertex for surface being shaded
in vec3 vert;

//vertex normal
in vec3 normal;

//directional vector of light
in vec3 eone;

//center point of area light
in vec3 lpos;

//radius of area light
in float rone;

//intensity?
uniform vec4 pointlight;




out vec4 color;

void main() {
    vec3 bright = vec3(0,0,0);
    float pi = 3.1415926535897932384626433832795;
    float incdegree = (2 * pi) / 360;


    vec3 z = vec3(0,0,1);

    //vcolor copied from pass.f.glsl
    vec3 colorrgb = vec3(115.0/255.0, 161.0/255.0, 227.0/255.0);

    for (float r = 0; r < rone; r+= 0.1){
        for(float theta = 0; theta < 2 * pi; theta += incdegree){
            vec3 pp = vec3(r * cos(theta),r * sin(theta),0);
            vec3 plane = cross(eone / length(eone), z / length(z));
            //Need to normalize even after the cross product of two unit vectors
            plane = normalize(plane);

            //rotation angle
            float rot = acos( dot(eone/length(eone),z / length(z)) / (
                length(eone)*length(z)) );

                //pp * vcos(rot) + (cross(plane, pp)*sin(rot)) + 
                //plane * (plane * pp)(1 - cos(rot))

            //rotation vector of pp + transform (lpos)
            vec3 Wpp = vec3(
                pp * cos(rot)  
                + normalize(cross(plane/length(plane), pp/ length(pp)))
                + plane * (dot(plane / length(plane),pp / length(pp))) * (1-cos(rot)) )
            + lpos;

            float tmp = pointlight[3] * dot(-lpos / length(lpos)), eone / length(eone)) 
                * (pi * (pow(rone,2)));
            //pointlight[3] = inherent light intensity 
            //* area * (-lpos * eone) to get light intensity

            Wpp *= tmp;
            colorrgb *= Wpp;

            //shade pixel
            //add shaded region to bright
            bright += colorrgb;

            color = (bright, 1);



            //Intensity of l = intensity of l knot * (-L * eone) * area (pi * radius * radius)

        }
    }
}