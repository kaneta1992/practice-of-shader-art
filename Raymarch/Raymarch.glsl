#define MAT_FLOOR 0.0

float sdPlane(vec3 p)
{
	return p.y;
}

vec2 opU(vec2 d1, vec2 d2)
{
	return (d1.x<d2.x) ? d1 : d2;
}

vec2 map(vec3 p)
{
    return vec2(sdPlane(p), MAT_FLOOR);
}

vec3 materialize(vec3 p, vec3 ray, float depth, vec2 mat)
{
    vec3 col = vec3(0.0);
    if (depth > 200.0) {
        col = vec3(0.0);
    } else if (mat.y == MAT_FLOOR) {
        float checker = mod(floor(p.x) + floor(p.z), 2.0);
        col = vec3(max(0.5, checker));
    }
    return col;
}

vec3 trace(vec3 p, vec3 ray)
{
    float t = 0.0;
    vec3 pos;
    vec2 mat;
    for (int i = 0; i < 60; i++) {
        pos = p + ray * t;
        mat = map(pos);
        if (mat.x < 0.001) {
        	break;
        }
        t += mat.x;
    }
    return materialize(pos, ray, t, mat);
}

mat3 camera(vec3 ro, vec3 ta, float cr )
{
	vec3 cw = normalize(ta - ro);
	vec3 cp = vec3(sin(cr), cos(cr),0.);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv = normalize( cross(cu,cw) );
    return mat3( cu, cv, cw );
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 p = (fragCoord.xy * 2.0 - iResolution.xy) / min(iResolution.x, iResolution.y);
    
    vec3 ro = vec3(0.0, 4.0, 0.0);
    vec3 ta = vec3(0.0, 40.0, 100.0);
    mat3 c = camera(ro, ta, 0.0);
    vec3 ray = c * normalize(vec3(p, 2.0));
    vec3 col = trace(ro, ray);

    fragColor = vec4(col,1.0);
}