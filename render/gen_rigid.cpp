#include<bits/stdc++.h>
using namespace std;

const int N=1e5+5,N2=105;

const float pi=acos(-1);
const float eps=1e-4;
const float g=10;
const float mu=0.001;
const float rho0=1e3;
const float k=40;
const float dt=2e-3;
// const float dt=1e-4;
const float h=0.04;
const float omega=0.05;
const float gamma=0.7;

const int M=1.3/h/2;

inline float sqr(float x){return x*x;}
inline float cubic(float x){return x*x*x;}
struct vec3{
	float x,y,z;
	inline float&operator[](int i){return i==0?x:i==1?y:z;}
};
inline float dot(vec3 a,vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
inline vec3 operator-(vec3 a,vec3 b){return vec3{a.x-b.x,a.y-b.y,a.z-b.z};}
inline vec3 operator+(vec3 a,vec3 b){return vec3{a.x+b.x,a.y+b.y,a.z+b.z};}
inline vec3 operator*(vec3 a,float b){return vec3{a.x*b,a.y*b,a.z*b};}
inline float len2(vec3 a){return dot(a,a);}
inline float len(vec3 a){return sqrt(len2(a));}
inline vec3 norm(vec3 a){return a*(1.0/len(a));}
inline vec3 cross(vec3 a,vec3 b){return vec3{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
inline float getW(vec3 r){
	float q=len(r)/h;
	if(q<=1){
		return 3/(2*pi*h*h*h)*(2.0/3-q*q+0.5*q*q*q);
	}else if(q<=2){
		return 3/(2*pi*h*h*h)*1.0/6*(2-q)*(2-q)*(2-q);
	}else return 0;
}
inline vec3 getWG(vec3 r){
	if(len(r)<eps)return vec3{0,0,0};
	float q=len(r)/h;
	vec3 res=r*(1.0/(len(r)*h));
	if(q<=1){
		return res*(3/(2*pi*h*h*h)*(-2*q+1.5*q*q));
	}else if(q<=2){
		return res*(3/(2*pi*h*h*h)*(-1.0/2*(2-q)*(2-q)));
	}else return vec3{0,0,0};
}

const vec3 BLUE={0.1, 0.6, 0.9};
const vec3 YELLOW={0.9, 0.5, 0.1};
const vec3 WHITE={1,1,1};
const vec3 colors[3]={BLUE,YELLOW,WHITE};

struct mat3{
	float a[3][3];
};
inline mat3 operator*(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=0;
		for(int k=0;k<3;++k){
			c.a[i][j]+=a.a[i][k]*b.a[k][j];
		}
	}
	return c;
}
inline vec3 operator*(mat3 a,vec3 b){
	vec3 c;
	for(int i=0;i<3;++i){
		c[i]=0;
		for(int j=0;j<3;++j){
			c[i]+=a.a[i][j]*b[j];
		}
	}
	return c;
}
inline mat3 operator+(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]+b.a[i][j];
	}
	return c;
}
inline mat3 operator*(mat3 a,float b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]*b;
	}
	return c;
}
inline mat3 operator-(mat3 a,mat3 b){
	mat3 c;
	for(int i=0;i<3;++i)for(int j=0;j<3;++j){
		c.a[i][j]=a.a[i][j]-b.a[i][j];
	}
	return c;
}
inline mat3 getinv(mat3 a){
	mat3 b;
	float det=a.a[0][0]*a.a[1][1]*a.a[2][2]+a.a[0][1]*a.a[1][2]*a.a[2][0]+a.a[0][2]*a.a[1][0]*a.a[2][1]
		-a.a[0][0]*a.a[1][2]*a.a[2][1]-a.a[0][1]*a.a[1][0]*a.a[2][2]-a.a[0][2]*a.a[1][1]*a.a[2][0];
	assert(det!=0);
	b.a[0][0]=(a.a[1][1]*a.a[2][2]-a.a[1][2]*a.a[2][1])/det;
	b.a[0][1]=(a.a[0][2]*a.a[2][1]-a.a[0][1]*a.a[2][2])/det;
	b.a[0][2]=(a.a[0][1]*a.a[1][2]-a.a[0][2]*a.a[1][1])/det;
	b.a[1][0]=(a.a[1][2]*a.a[2][0]-a.a[1][0]*a.a[2][2])/det;
	b.a[1][1]=(a.a[0][0]*a.a[2][2]-a.a[0][2]*a.a[2][0])/det;
	b.a[1][2]=(a.a[0][2]*a.a[1][0]-a.a[0][0]*a.a[1][2])/det;
	b.a[2][0]=(a.a[1][0]*a.a[2][1]-a.a[1][1]*a.a[2][0])/det;
	b.a[2][1]=(a.a[0][1]*a.a[2][0]-a.a[0][0]*a.a[2][1])/det;
	b.a[2][2]=(a.a[0][0]*a.a[1][1]-a.a[0][1]*a.a[1][0])/det;
	return b;
}
inline mat3 a2m(vec3 a){
	mat3 b;
	b.a[0][0]=0;
	b.a[0][1]=-a.z;
	b.a[0][2]=a.y;
	b.a[1][0]=a.z;
	b.a[1][1]=0;
	b.a[1][2]=-a.x;
	b.a[2][0]=-a.y;
	b.a[2][1]=a.x;
	b.a[2][2]=0;
	return b;
}
const mat3 I3=mat3{1,0,0,0,1,0,0,0,1};

vec3 center[N2],initial_center[N2];
mat3 trans[N2];
bool bo[N2];

int main(){
    for(int i=0;i<600;++i){
        char c[99];
        sprintf(c,"raw/%d_rigid.txt",i);
        FILE*fi=fopen(c,"r");
        int rcnt;
        fscanf(fi,"%d",&rcnt);
        for(int j=1;j<=rcnt;++j){
            fscanf(fi,"%f%f%f",&center[j].x,&center[j].y,&center[j].z);
			if(!bo[j]){
				bo[j]=1;
				initial_center[j]=center[j];
			}
            for(int k=0;k<3;++k){
                fscanf(fi,"%f%f%f",&trans[j].a[k][0],&trans[j].a[k][1],&trans[j].a[k][2]);
            }
        }
        fclose(fi);
        FILE*fo[2];
        sprintf(c,"tmp3/rigid_yellow_%d.obj",i);
        fo[0]=fopen(c,"w");
        sprintf(c,"tmp3/rigid_green_%d.obj",i);
        fo[1]=fopen(c,"w");
		int id[30];
		for(int i=1;i<=3;++i)id[i]=0;
		id[4]=id[6]=id[8]=id[9]=id[11]=id[13]=id[14]=id[16]=id[18]=1;
		id[5]=id[7]=id[10]=id[12]=id[15]=id[17]=0;
		for(int i=19;i<=27;++i)id[i]=i&1;
		assert(rcnt<=27);
		int cnt[2]={0,0};
        for(int j=1;j<=rcnt;++j){
            sprintf(c,"tmp2/output%d.obj",j);
            fi=fopen(c,"r");
			int oldc=cnt[id[j]];
			for(;fscanf(fi,"%s",c)!=EOF;){
				if(*c=='v'){
					vec3 v;
					fscanf(fi,"%f%f%f",&v.x,&v.y,&v.z);
					v=trans[j]*(v-initial_center[j])+center[j];
					fprintf(fo[id[j]],"v %f %f %f\n",v.x,v.z,v.y);
					++cnt[id[j]];
				}else{
					int a,b,c;
					fscanf(fi,"%d%d%d",&a,&b,&c);
					fprintf(fo[id[j]],"f %d %d %d\n",oldc+a,oldc+b,oldc+c);
				}
			}
			fclose(fi);
        }
		for(int i=0;i<2;++i){
			if(cnt[i]==0){
				fprintf(fo[i],"v 9 10 10\nv 10 11 10\nv 10 10 11\nv 10 11 11\nf 1 2 3\nf 1 3 4\nf 1 4 2\nf 4 3 2\n");
			}
			fclose(fo[i]);
		}
    }
}