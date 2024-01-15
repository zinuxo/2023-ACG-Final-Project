#include<bits/stdc++.h>
using namespace std;

const int N=2e5+5;

const float pi=acos(-1);
const float eps=1e-4;
const float g=10;
const float mu=0.0005;
const float rho0=1e3;
const float k=40;
const float dt=5e-3;
const float h=0.04;

const int M=1.3/h/2;

inline float sqr(float x){return x*x;}
inline float cubic(float x){return x*x*x;}
struct vec3{float x,y,z;};
float dot(vec3 a,vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
vec3 operator-(vec3 a,vec3 b){return vec3{a.x-b.x,a.y-b.y,a.z-b.z};}
vec3 operator+(vec3 a,vec3 b){return vec3{a.x+b.x,a.y+b.y,a.z+b.z};}
vec3 operator*(vec3 a,float b){return vec3{a.x*b,a.y*b,a.z*b};}
float len(vec3 a){return sqrt(dot(a,a));}
vec3 norm(vec3 a){return a*(1.0/len(a));}

struct particle{
	vec3 x,v,vv,a;
	float m,rho;
	int type; // 0 dynamic, 1 static
}a[N];
int n;
void add(int type,int nx,int ny,int nz,float xs,float xl,float ys,float yl,float zs,float zl){
	float r=h/2;
	for(int i=1;i<=nx;++i)
		for(int j=1;j<=ny;++j)
			for(int k=1;k<=nz;++k){
				a[++n]=particle{
					vec3{xs+xl*i/(nx+1),ys+yl*j/(ny+1),zs+zl*k/(nz+1)},
					vec3{0,0,-g/2},
					vec3{0,0,-g/2},
					vec3{0,0,-g},
					pi*4/3*r*r*r*rho0,rho0,type
				};
	}
}
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

vec3 lb,rb;
struct hmp{
	inline void getB(){
		lb=rb=a[1].x;
		for(int i=2;i<=n;++i){
			lb.x=min(lb.x,a[i].x.x);
			rb.x=max(rb.x,a[i].x.x);
			lb.y=min(lb.y,a[i].x.y);
			rb.y=max(rb.y,a[i].x.y);
			lb.z=min(lb.z,a[i].x.z);
			rb.z=max(rb.z,a[i].x.z);
		}
	}
	int be[M*M*M];
	particle b[N];
	inline void rebuild(){
		int i;
		memset(be,0,sizeof be);
		for(i=1;i<=n;++i){
			int x=floor((a[i].x.x-lb.x+eps)/h/2);
			int y=floor((a[i].x.y-lb.y+eps)/h/2);
			int z=floor((a[i].x.z-lb.z+eps)/h/2);
			if(z>=M || z<0){
				cerr<<i<<' '<<a[i].x.z<<' '<<lb.z<<' '<<h<<endl;
			}
			assert(0<=x && x<M);
			assert(0<=y && y<M);
			assert(0<=z && z<M);
			be[x*M*M+y*M+z]++;
		}
		int max=0;
		for(i=1;i<M*M*M;++i){
			max=max>be[i]?max:be[i];
			be[i]+=be[i-1];
		}
		// cerr<<max<<endl;
		for(i=n;i>=1;--i){
			int x=floor((a[i].x.x-lb.x)/h/2);
			int y=floor((a[i].x.y-lb.y)/h/2);
			int z=floor((a[i].x.z-lb.z)/h/2);
			b[be[x*M*M+y*M+z]--]=a[i];
		}		
	}
}H;

void step(){
	for(int i=1;i<=n;++i)if(a[i].type==0){
		a[i].x=a[i].x+a[i].vv*dt;
		a[i].v=a[i].vv+a[i].a*(dt/2);
	}

	for(int i=1;i<=n;++i)if(a[i].type==0){
		if(a[i].x.x<0+h){
			a[i].x.x=0+h;
			a[i].v.x*=-0.5;
			a[i].vv.x*=-0.5;
			a[i].a.x*=-0.5;
		}
		if(a[i].x.x>1-h){
			a[i].x.x=1-h;
			a[i].v.x*=-0.5;
			a[i].vv.x*=-0.5;
			a[i].a.x*=-0.5;
		}
		if(a[i].x.y<0+h){
			a[i].x.y=0+h;
			a[i].v.y*=-0.5;
			a[i].vv.y*=-0.5;
			a[i].a.y*=-0.5;
		}
		if(a[i].x.y>1-h){
			a[i].x.y=1-h;
			a[i].v.y*=-0.5;
			a[i].vv.y*=-0.5;
			a[i].a.y*=-0.5;
		}
		if(a[i].x.z<0+h){
			a[i].x.z=0+h;
			a[i].v.z*=-0.5;
			a[i].vv.z*=-0.5;
			a[i].a.z*=-0.5;
		}
		if(a[i].x.z>1-h){
			a[i].x.z=1-h;
			a[i].v.z*=-0.5;
			a[i].vv.z*=-0.5;
			a[i].a.z*=-0.5;
		}
	}

	H.rebuild();
	for(int i=1;i<=n;++i)if(a[i].type==0){
		float new_rho=0;
		int x=floor((a[i].x.x-lb.x)/h/2);
		int y=floor((a[i].x.y-lb.y)/h/2);
		int z=floor((a[i].x.z-lb.z)/h/2);

		for(int dx=-1;dx<=1;++dx)
			for(int dy=-1;dy<=1;++dy)
				for(int dz=-1;dz<=1;++dz){
					int xx=x+dx,yy=y+dy,zz=z+dz;
					if(xx<0 || xx>=M || yy<0 || yy>=M || zz<0 || zz>=M)continue;
					int id=xx*M*M+yy*M+zz;
					int start=H.be[id]+1,end=id==M*M*M-1?n:H.be[id+1]-1;
					for(int j=start;j<=end;++j){
						if(len(a[i].x-H.b[j].x)<=2*h){
							new_rho+=H.b[j].m*getW(a[i].x-H.b[j].x);
						}
					}
				}
				
		new_rho=max(new_rho,rho0);
		a[i].rho=new_rho;
	}
	H.rebuild();

	int tot=0;

	for(int i=1;i<=n;++i)if(a[i].type==0){
		vec3 f=vec3{0,0,-g*a[i].m};
		int x=floor((a[i].x.x-lb.x)/h/2);
		int y=floor((a[i].x.y-lb.y)/h/2);
		int z=floor((a[i].x.z-lb.z)/h/2);

		for(int dx=-1;dx<=1;++dx)
			for(int dy=-1;dy<=1;++dy)
				for(int dz=-1;dz<=1;++dz){
					int xx=x+dx,yy=y+dy,zz=z+dz;
					if(xx<0 || xx>=M || yy<0 || yy>=M || zz<0 || zz>=M)continue;
					int id=xx*M*M+yy*M+zz;
					int start=H.be[id]+1,end=id==M*M*M-1?n:H.be[id+1]-1;
					for(int j=start;j<=end;++j){
++tot;
						particle aj=H.b[j];
						if(len(a[i].x-aj.x)<=2*h && len(a[i].x-aj.x)>eps){
							// f=f+getWG(a[i].x-aj.x)*
							// 	(mu*10*(aj.m/aj.rho)*dot(a[i].v-aj.v,a[i].x-aj.x)
							// 	/sqr(len(a[i].x-aj.x)));
							f=f-(a[i].v-aj.v)*(mu*(aj.m/aj.rho)*2*len(getWG(a[i].x-aj.x))/len(a[i].x-aj.x));
							float pi=k*(a[i].rho-rho0);
							float pj=k*(aj.rho-rho0);
							f=f-getWG(a[i].x-aj.x)*aj.m*(pi/sqr(a[i].rho)+pj/sqr(aj.rho));
if(isnan(f.z))
++i,--i;
if(isinf(f.z) || isinf(f.x))
++i,--i;
if(fabs(f.z)>3000)
++i,--i,getWG(a[i].x-aj.x);
if((getWG(a[i].x-aj.x)*aj.m*(pi/sqr(a[i].rho)+pj/sqr(aj.rho))).x>eps)
++i,--i;
						}
					}
				}
		// vec3 old_a=a[i].a;
		// vec3 new_a=f*(1.0/a[i].m);
		// a[i].x=a[i].x+a[i].v*dt+old_a*(dt*dt/2);
		// a[i].v=a[i].v+(old_a+new_a)*(dt/2);
		// a[i].a=new_a;
		vec3 new_a=f*(1.0/a[i].m);
		a[i].a=new_a;
		a[i].vv=a[i].vv+a[i].a*dt;
	}
// cerr<<tot<<endl;
}

char obuf[1<<20],*oh;
inline void write(float x){
	if(x<0){
		putchar('-');
		x=-x;
	}
	int t=floor(x);
	assert(t<10);
	*oh++=t+'0';
	*oh++='.';
	x-=t;
	for(int i=0;i<5;++i){
		x*=10;
		t=floor(x);
		*oh++=t+'0';
		x-=t;
	}
	*oh++=' ';
}
inline void out(vec3 P){
	write(P.x);
	write(P.z);
	write(P.y);
	*oh++=' ';
	write(0.1);
	write(0.6);
	write(0.9);
	*oh++='\n';
}

int main(){
	int nn=10;
	add(0,nn,nn,nn,0.1,0.5,0.1,0.5,0.25,0.5);	
	
	float rr=h/2;
	add(1,40,40,1,-rr,1+rr*2,-rr,1+rr*2,0,0);
	add(1,40,40,1,-rr,1+rr*2,-rr,1+rr*2,-rr/2,0);
	add(1,40,40,1,-rr,1+rr*2,-rr,1+rr*2,1,0);
	add(1,40,40,1,-rr,1+rr*2,-rr,1+rr*2,1+rr/2,0);
	add(1,40,1,40,-rr,1+rr*2,0,0,-rr,1+rr*2);
	add(1,40,1,40,-rr,1+rr*2,-rr/2,0,-rr,1+rr*2);
	add(1,40,1,40,-rr,1+rr*2,1,0,-rr,1+rr*2);
	add(1,40,1,40,-rr,1+rr*2,1+rr/2,0,-rr,1+rr*2);
	add(1,1,40,40,0,0,-rr,1+rr*2,-rr,1+rr*2);
	add(1,1,40,40,-rr/2,0,-rr,1+rr*2,-rr,1+rr*2);
	add(1,1,40,40,1,0,-rr,1+rr*2,-rr,1+rr*2);
	add(1,1,40,40,1+rr/2,0,-rr,1+rr*2,-rr,1+rr*2);

	// add(2,6,6,6,0.8,0.1,0.8,0.1,0,0.1);

	H.getB();
	lb=vec3{-rr,-rr,-rr};
	rb=vec3{1,1,1};

	float tt=0.01;
	float rt=tt;
	for(int fr=0;fr<100;){
		step();
		rt-=dt;
		if(rt<=0){
			rt=tt;
			char c[99];
			sprintf(c,"raw/%d.txt",fr++);
			FILE*f=fopen(c,"w");
			oh=obuf;
			for(int i=1;i<=n;++i)if(a[i].type!=1){
				out(a[i].x);
			}
			for(int i=0;i<8;++i){
				// out(vec3{float(i>>2&1),float(i>>1&1),float(i&1)});
			}
			fwrite(obuf,1,oh-obuf,f);
			fclose(f);
			fprintf(stderr,"frame:%d time:%d\n",fr,clock());
		}
	}
}