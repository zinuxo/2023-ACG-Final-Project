#include"scene.hpp"

int main(){
	int nn=8;
	// add(0,nn,nn,nn,0.1,0.5,0.1,0.5,0.1,0.5);
	// add(0,nn*2,nn*2,nn,0.1,0.9,0.1,0.9,0.1,0.5);
	// add(0,nn*2,nn*2,nn,0.1,0.8,0.1,0.8,0.3,0.4);
	// add(0,nn*2,nn*2,nn/2,0.1,0.8,0.1,0.8,0.1,0.2);
	// add(0,nn,nn,nn,0.1,0.5,0.1,0.5,0.75,0.2);
	// add(0,16,16,16,0.1,0.8,0.1,0.8,0.1,0.8);
	// add(0,16,16,16,2.1,0.8,0.1,0.8,0.1,0.8);
	// add(0,48,32,16,0.3,2.4,0.2,1.6,1.1,0.8);
	// add(0,48,32,16,0.2,2.4,0.2,1.6,0.9,0.8);
	add(0,60,40,10,0.2,2.4,0.2,1.6,0.1,0.4);
	
	float rr=h/2;
	// add(1,20,20,1,-rr,1+rr*2,-rr,1+rr*2,0,0);
	// add(1,20,20,1,-rr,1+rr*2,-rr,1+rr*2,-rr/2,0);
	// add(1,20,20,1,-rr,1+rr*2,-rr,1+rr*2,1,0);
	// add(1,20,20,1,-rr,1+rr*2,-rr,1+rr*2,1+rr/2,0);
	// add(1,20,1,20,-rr,1+rr*2,0,0,-rr,1+rr*2);
	// add(1,20,1,20,-rr,1+rr*2,-rr/2,0,-rr,1+rr*2);
	// add(1,20,1,20,-rr,1+rr*2,1,0,-rr,1+rr*2);
	// add(1,20,1,20,-rr,1+rr*2,1+rr/2,0,-rr,1+rr*2);
	// add(1,1,20,20,0,0,-rr,1+rr*2,-rr,1+rr*2);
	// add(1,1,20,20,-rr/2,0,-rr,1+rr*2,-rr,1+rr*2);
	// add(1,1,20,20,1,0,-rr,1+rr*2,-rr,1+rr*2);
	// add(1,1,20,20,1+rr/2,0,-rr,1+rr*2,-rr,1+rr*2);
	int oldn=A.n;
	for(int i=0;i<3;++i){
		int nnn[3]={60,40,40};
		float ss[3]={-rr,-rr,-rr};
		float le[3]={3+rr*2,2+rr*2,2+rr*2};
		nnn[i]=1;
		le[i]=0;
		add(1,nnn[0],nnn[1],nnn[2],ss[0],le[0],ss[1],le[1],ss[2],le[2]);
		ss[i]-=rr/2;
		add(1,nnn[0],nnn[1],nnn[2],ss[0],le[0],ss[1],le[1],ss[2],le[2]);
		ss[i]=i==0?3:2;
		add(1,nnn[0],nnn[1],nnn[2],ss[0],le[0],ss[1],le[1],ss[2],le[2]);
		ss[i]+=rr/2;
		add(1,nnn[0],nnn[1],nnn[2],ss[0],le[0],ss[1],le[1],ss[2],le[2]);
	}
	for(int i=oldn+1;i<=A.n;++i){
		float sum=0;
		A.rho[i]=rho0*3;
		for(int j=oldn+1;j<=A.n;++j)if(len(A.x[i]-A.x[j])<=h*2)
			sum+=getW(A.x[i]-A.x[j]);
		A.mass[i]=gamma/sum*A.rho[i];
	}
	auto getv=[&](int i){
		return vec3{cos(i*2*pi/6),sin(i*2*pi/6),0};
	};

	float hh=0.8;
	vec3 o1=vec3{1.5,1.6,hh},o2=vec3{1.5,1.6,hh}+getv(4)*1.32,o3=vec3{1.5,1.6,hh}+getv(5)*1.32;
	add2(o1,0.25,{getv(4),getv(5)});
	add2(o2,0.25,{getv(0),getv(1)});
	add2(o3,0.25,{getv(2),getv(3)});
	auto add6=[&](vec3 o1,vec3 o2){
		vec3 p=norm(o2-o1);
		vec3 q=vec3{-p.y,p.x,p.z};
		vector<int>tmp={34,35,36,37,38};
		for(int i=0;i<5;++i)if(i&1){
			add3(o1+p*(0.42+i*0.12),p,vec3{0,0,1},1,tmp);
		}else{
			add3(o1+p*(0.42+i*0.12),p,q,3,tmp);
		}
	};
	add6(o1,o2);
	add6(o2,o3);
	add6(o3,o1);


	// for(int i=oldn+1;i<=A.n;++i){
	// 	A.info[i]=1;
	// 	A.v[i].x=114514;
	// }
	// --rcnt;


	for(float x:{0.2,0.45,0.7})for(float y:{0.2,0.45,0.7}){		
		// add(2,4,4,4,x,0.15,y,0.15,0.55,0.15);
	}
	// add(2,4,4,4,2,0.15,0.4,0.15,0.55,0.15);
	// add(2,4,4,4,0.4,0.15,0.4,0.15,0.55,0.15);
	// add(2,4,4,4,0.6,0.15,0.6,0.15,0.55,0.15);
	// add(2,4,4,4,0.8,0.15,0.8,0.15,0.55,0.15);
	// add(2,4,4,4,0.5,0.15,0.5,0.15,0.75,0.15);
	// add(2,4,4,4,0.45,0.15,0.45,0.15,0.75,0.15);
	// add(2,4,4,4,0.4,0.15,0.4,0.15,0.1,0.15);
	// add(2,4,4,4,0.6,0.15,0.6,0.15,0.1,0.15);
	// add(2,4,4,4,0.8,0.15,0.8,0.15,0.1,0.15);

	auto add5=[&](vec3 o,float r,function<bool(vec3,float)>ck,int color=1){
		++rcnt;
		int oldn=A.n;

		float r0=0.03;
		int t=ceil(r/r0)*2;
		for(int i=-t;i<=t;++i)for(int j=-t;j<=t;++j)for(int k=-t;k<=t;++k){
			vec3 p=vec3{i,j,k}*r0;
			if(!ck(p,r))continue;
			++A.n;
			A.x[A.n]=o+p;
			A.v[A.n]=vec3{0,0,0};
			A.rho[A.n]=rho0;
			A.info[A.n]=2+(rcnt<<8)+(color<<4);
		}

		ri[rcnt].mass=0;
		ri[rcnt].I=I3*0;
		ri[rcnt].c=o;
		ri[rcnt].v=vec3{0,0,0};
		ri[rcnt].w=vec3{0,0,0};
		for(int i=oldn+1;i<=A.n;++i){
			vec3 p=A.x[i]-o;
			mat3 dI=mat3{
				sqr(p.y)+sqr(p.z),-p.x*p.y,-p.x*p.z,
				-p.x*p.y,sqr(p.x)+sqr(p.z),-p.y*p.z,
				-p.x*p.z,-p.y*p.z,sqr(p.x)+sqr(p.y)
			};
			ri[rcnt].I=ri[rcnt].I+dI;
			ri[rcnt].mass=ri[rcnt].mass+rho0/5*r0*r0*r0;
		}
		ri[rcnt].I=ri[rcnt].I*ri[rcnt].mass;
		ri[rcnt].Iinv=getinv(ri[rcnt].I);
		ri[rcnt].trans=I3;
		for(int i=oldn+1;i<=A.n;++i){
			float sum=0;
			for(int j=oldn+1;j<=A.n;++j)if(len(A.x[i]-A.x[j])<=h*2)
				sum+=getW(A.x[i]-A.x[j]);
			A.mass[i]=gamma/sum*A.rho[i];
		}
	};
	auto ck1=[&](vec3 p,float d){
		float y=p.z;
		float x=sqrt(sqr(p.x)+sqr(p.y));
		return pow(fabs(x),2.0/3)+pow(fabs(y),2.0/3)<=pow(d,2.0/3);
	};
	auto ck2=[&](vec3 p,float d){
		// float y=p.z;
		// float x=sqrt(sqr(p.x)+sqr(p.y));
		// float r=sqrt(sqr(x)+sqr(y));
		// float theta=atan2(y,x);
		// return r<=d*(1-sin(theta));
		return 3*sqr(p.x)+3*sqr(p.y)+sqr(p.z)<=sqr(d);
	};
	auto ck3=[&](vec3 p,float d){
		return fabs(p.x)+fabs(p.y)+fabs(p.z)<=d;
	};
	auto work=[&](vec3 o,int color){
		float r=0.07;
		add5(o,r,ck3,color);
		char c[99];
		sprintf(c,"tmp2/output%d.obj",rcnt);
		freopen(c,"w",stdout);
		vec3 x;
		x=o+vec3{0,0,r};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		x=o+vec3{0,0,-r};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		x=o+vec3{0,r,0};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		x=o+vec3{0,-r,0};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		x=o+vec3{r,0,0};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		x=o+vec3{-r,0,0};
		printf("v %.6f %.6f %.6f\n",x.x,x.y,x.z);
		for(int i:{1,2})for(int j:{3,4})for(int k:{5,6})printf("f %d %d %d\n",i,j,k);
	};
	for(int i=0;i<3;++i){
		// function<bool(vec3,float)>cks[3]={ck1,ck2,ck3};
		for(int j=0;j<4;++j){
			// work(vec3{0.5+j*0.5,0.5+i*0.5,1.3});
		}
	}
	cerr<<"gua\n";
	float tt=0.005;
	float rt=0;
	for(int fr=0;fr<600;){
		rt-=dt;
		if(rt<=0){
			if(fr%4==0 && fr<1000 && fr>100){
				int oldn=A.n;
				add(0,8,8,2,1.3,0.4,0.9,0.4,1.6,0.15);
				for(int i=oldn+1;i<=A.n;++i)A.v[i]=vec3{0,0,-2},time_p[i]=fr;
			}
			if(fr%40==15 && fr<480 && fr>100){
				int oldn=A.n;
				work(vec3{1.5,0.4,1.6},fr/20%2==1?1:3);
				for(int i=oldn+1;i<=A.n;++i)time_p[i]=fr;
				time_r[rcnt]=fr;
			}
			++fr;
		}
	}
    freopen("scene.txt","w",stdout);
    printf("%d\n",A.n);
    for(int i=1;i<=A.n;++i){		
		printf("%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d %d\n",
			A.x[i].x,A.x[i].y,A.x[i].z,
			A.v[i].x,A.v[i].y,A.v[i].z,
			A.mass[i],A.rho[i],A.info[i],time_p[i]
		);
	}
    printf("%d\n",rcnt);
    for(int i=1;i<=rcnt;++i){
        ri[i].output();
        printf("%d\n",time_r[i]);
    }
	printf("600\n0.005\n");
}