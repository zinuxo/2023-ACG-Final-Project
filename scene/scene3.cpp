#include"scene.hpp"

int main(){
	add(0,55,35,45,0.1,2,0.3,1.5,0.1,1.8);
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
	printf("50\n0.005\n");
}