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

// iq's sky https://www.shadertoy.com/view/MdX3Rr
vec3 sunDir = normalize(vec3(.2, .25, .5));
vec3 skyColor(vec3 rd)
{
    float sundot = clamp(dot(rd,sunDir),0.0,1.0);
    // sky		
    vec3 col = vec3(0.2,0.5,0.85)*1.1 - rd.y*rd.y*0.5;
    col = mix( col, 0.85*vec3(0.7,0.75,0.85), pow( 1.0-max(rd.y,0.0), 4.0 ) );
    // sun
    col += 0.25*vec3(1.0,0.7,0.4)*pow( sundot,5.0 );
    col += 0.25*vec3(1.0,0.8,0.6)*pow( sundot,64.0 );
    col += 0.2*vec3(1.0,0.8,0.6)*pow( sundot,512.0 );
    // horizon
    col = mix( col, 0.68*vec3(0.4,0.65,1.0), pow( 1.0-max(rd.y,0.0), 16.0 ) );
    return col;
}

vec3 materialize(vec3 p, vec3 ray, float depth, vec2 mat)
{
    vec3 col = vec3(0.0);
    vec3 sky = skyColor(ray);
    if (depth > 2000.0) {
        col = sky;
    } else if (mat.y == MAT_FLOOR) {
        float checker = mod(floor(p.x) + floor(p.z), 2.0);
        col = vec3(max(0.5, checker)) * vec3(0.4,0.6,0.85) * 1.1;
    }
    float fo = 1.0-exp(-pow(0.002*depth,1.5) );
    vec3 fco = 0.65*vec3(0.4,0.65,1.0);
    col = mix( col, sky, fo );
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
    
    // flare
    float sundot = clamp(dot(ray,sunDir),0.0,1.0);
    col += 0.3*vec3(1.0,0.7,0.3)*pow( sundot, 8.0 );
    col = pow(col, vec3(1.0/2.2));

    fragColor = vec4(col,1.0);
}