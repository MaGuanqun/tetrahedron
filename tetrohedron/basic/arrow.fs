#version 460 core
out vec4 FragColor;

//in vec2 TexCoord;
in vec3 Normal;
in vec3 FragPos;
//in vec3 lighting;

uniform vec3 viewPos;
uniform vec3 lightPos;
uniform int type;

void main()
{ 
    if(type==0)
        FragColor =vec4 (1.0, 0.1, 0.1, 1.0);
    else if(type==1)
        FragColor =vec4(0.1,1.0,0.1,1.0);
    else
        FragColor=vec4(0.1,0.1,1.0,1.0);    
    vec3 ambient=vec3(0.4,0.4,0.4);

    vec3 norm=normalize(Normal);
    vec3 lightDir=normalize(lightPos-FragPos);
    float diff=max(dot(norm,lightDir),0.0);
    vec3 diffuse=diff*vec3(0.8,0.8,0.8);

    vec3 viewDir=normalize(viewPos-FragPos);
    vec3 reflectDir=reflect(-lightDir,norm);
    float spec=pow(max(dot(viewDir,reflectDir),0.0),32);
    vec3 specular=spec*vec3(0.95,0.95,0.95);

    //float distance=sqrt(dot(light.position-FragPos,light.position-FragPos));
    //float attenuation=1.0f/(light.constant+light.linear*distance+light.quadratic*distance*distance);

   // diffuse*=attenuation;
   // specular*=attenuation;

    vec3 result=ambient+diffuse+specular;
 
    FragColor*=vec4(result,1.0f);
    
       

}