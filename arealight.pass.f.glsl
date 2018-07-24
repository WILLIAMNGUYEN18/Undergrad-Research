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


out vec4 color;

void main() {
    float bright = 0;
    incdegree = (2 * 3.1415926535897932384626433832795) / 360;

    vec3 z = vec3(0,0,1);

    for (float r = 0; r < rone; r+= 0.1){
        for(float theta = 0; theta < incdegree * 360; theta += incdegree){
            vec3 pp = vec(r * cos(theta),r * sin(theta),0);
            vec3 plane = cross(eone / length(eone), z / length(z));
            //Need to normalize even after the cross product of two unit vectors
            plane = normalize(plane);

            //rotation angle
            float rot = acos( dot(eone/length(eone),z / length(z)) / (
                length(eone)*length(z)) );

                //pp * vcos(rot) + (cross(plane, pp)*sin(rot)) + 
                //plane * (plane * pp)(1 - cos(rot))

            //
            vec3 Wpp = vec(pp * cos(rot) + 
            normalize(cross(plane/length(plane), pp/ length(pp))) +
            plane * (dot(plane / length(plane),pp / length(pp))) * (1-cos(rot)) ) +
            lpos; 
            //Intensity of l = intensity of l knot * (-L * eone) * area (pi * radius * radius)

        }
    }
}