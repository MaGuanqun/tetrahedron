#version 430 core

out vec4 FragColor;

in VS_OUT{
	vec3 FragPos;
	vec3 Normal;
} fs_in;

struct Light{
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;
};

uniform vec3 viewPos;
uniform vec3 lightPos;

uniform samplerCube depthMap;
uniform float far_plane;
uniform bool lightShadowOn;
uniform bool lightIsChosen;
uniform Light light;

float PointShadowCalculation(vec3 FragPos)
{
    vec3 fragToLight = FragPos - lightPos;
    // get vector between the fragment position     
    // now get current linear depth as the length between the fragment and light position
    float currentDepth = length(fragToLight);
    float bias = 0.05f; // we use a much larger bias since depth is now in [near_plane, far_plane] range
    float shadow = 0.0f;      
	int samples=20;
	float offset=0.01f;
    float viewDistance = length(viewPos - FragPos);
    float diskRadius = (1.0f + viewDistance / far_plane) / 200.0f; 

    vec3 sampleOffsetDirections[20] = vec3[]
    (
       vec3( 1,  1,  1), vec3( 1, -1,  1), vec3(-1, -1,  1), vec3(-1,  1,  1), 
       vec3( 1,  1, -1), vec3( 1, -1, -1), vec3(-1, -1, -1), vec3(-1,  1, -1),
       vec3( 1,  1,  0), vec3( 1, -1,  0), vec3(-1, -1,  0), vec3(-1,  1,  0),
       vec3( 1,  0,  1), vec3(-1,  0,  1), vec3( 1,  0, -1), vec3(-1,  0, -1),
       vec3( 0,  1,  1), vec3( 0, -1,  1), vec3( 0, -1, -1), vec3( 0,  1, -1)
    );   
	for(int i = 0; i < samples; ++i)
    {
		float closeDepth=texture(depthMap,fragToLight+sampleOffsetDirections[i]*diskRadius).r;
	    closeDepth*=far_plane;
		if(currentDepth-bias>closeDepth)
	    {
			shadow+=1.0f;
		}		
    }
    shadow/=float(samples);
    return shadow;
}


void main()
{
	FragColor=vec4(0.02,0.02,0.6,1.0);	
    if(lightIsChosen){
        float constant = 0.06f;
        float linear = 0.09f;
        float quadratic = 0.05f;   
	    vec3 norm=normalize(fs_in.Normal);
	    vec3 lightDir=normalize(lightPos-fs_in.FragPos);

        float diff=max(dot(norm,lightDir),0.0);
        vec3 diffuse=diff*light.diffuse;

	    vec3 viewDir=normalize(viewPos-fs_in.FragPos);	
	    vec3 reflectDir=reflect(-lightDir,norm);
	    float spec=pow(max(dot(viewDir,reflectDir),0.0),32);
        vec3 specular=spec*light.specular;

        //float distance=sqrt(dot(lightPos-fs_in.FragPos,lightPos-fs_in.FragPos));
       // float attenuation =1.0f/(constant+linear*distance+quadratic*distance*distance);

       // diffuse*=attenuation;
       // specular*=attenuation;
       vec3 result;
       float shadow;
       if(lightShadowOn){
       shadow=PointShadowCalculation(fs_in.FragPos);
        result=light.ambient+(1.0-shadow)*(diffuse+specular);
       }
       else{
        result=light.ambient+diffuse+specular;
        }
       FragColor*=vec4(result,1.0f);
   }
}