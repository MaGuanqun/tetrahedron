#version 460 core

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

struct Material{
    vec3 Kd;
    vec3 Ka;
    vec3 Ks;
    float Ns;
};

uniform vec3 viewPos;
uniform vec3 lightPos;

uniform samplerCube depthMap;
uniform float far_plane;
uniform bool lightShadowOn;
uniform bool lightIsChosen;
uniform Light light;
uniform Material material;
uniform float transparence;


float PointShadowCalculation(vec3 FragPos)
{
    vec3 fragToLight = FragPos - lightPos;
    float currentDepth = length(fragToLight);
    float bias =0.01*far_plane; // we use a much larger bias since depth is now in [near_plane, far_plane] range
    float shadow = 0.0;      
    float samples = 4.0;
    float offset  = 0.001*far_plane;
    for(float x = -offset; x < offset; x += offset / (samples * 0.5))
    {
        for(float y = -offset; y < offset; y += offset / (samples * 0.5))
        {
            for(float z = -offset; z < offset; z += offset / (samples * 0.5))
            {
                float closestDepth = texture(depthMap, fragToLight + vec3(x, y, z)).r; 
                closestDepth *= far_plane;   // undo mapping [0;1]
                if(currentDepth - bias > closestDepth)
                    shadow += 1.0;
            }
        }
    }
    shadow/=float(samples*samples*samples);
   return shadow;
}



void main()
{
    FragColor=vec4(material.Kd,1.0);	
    if(lightIsChosen){
        //float constant = 0.06f;
        // float linear = 0.09f;
        //float quadratic = 0.05f;   


	    vec3 norm=normalize(fs_in.Normal);
	    vec3 lightDir=normalize(lightPos-fs_in.FragPos);
        float diff=abs(dot(norm,lightDir));
        vec3 diffuse=diff*light.diffuse*material.Kd;

	    vec3 viewDir=normalize(viewPos-fs_in.FragPos);	
	    vec3 reflectDir=reflect(-lightDir,norm);
	    float spec=pow(abs(dot(viewDir,reflectDir)),material.Ns);
        vec3 specular=spec*light.specular*material.Ks;
        //float distance=sqrt(dot(lightPos-fs_in.FragPos,lightPos-fs_in.FragPos));
        // float attenuation =1.0f/(constant+linear*distance+quadratic*distance*distance);
        // diffuse*=attenuation;
        // specular*=attenuation;
        vec3 result;
        float shadow;
        if(lightShadowOn){
        shadow=PointShadowCalculation(fs_in.FragPos);
        //float shadow=ShadowCalculation(norm);
        result=light.ambient*material.Ka+(1.0-shadow)*(diffuse+specular);
        }
        else{
        result=light.ambient*material.Ka+diffuse+specular;
        }
        FragColor=vec4(result,transparence);
        }

}