#include "stdio.h"
#include "math.h"

// this is the function that does the work. declare here, define later
double richardsonDistance(double* X, double* Y, int N, double dist, double* polygon_pts, double* polyArea);
void linInterp(double x1, double x2,double y1,double y2,double *xstar,double* ystar,double dist);

int main(void) {
// define file pointers to open files with
FILE* boundary_fp;
FILE* distance_fp;
FILE* params_fp;
FILE* richdist_fp;
FILE* polypts_fp;
// define integers for lengths of the boundary and number of distance
// points to use, as well as the file I/O return code rc and index var i.
int b_len,d_len,rc,i,j;
// read the parameter file to get the lengths of the boundary and distance arrays
params_fp = fopen("params.dat","r");
rc=fscanf(params_fp,"%d",&b_len);
rc=fscanf(params_fp,"%d",&d_len);
fclose(params_fp);
// initialize arrays now that we have their sizes
double boundary_x[b_len];
double boundary_y[b_len];
double distance[d_len];
double richdist[d_len];
double polygon_pts[2*b_len];
double polyarea=0;
// read in the boundary data
boundary_fp = fopen("boundary.dat","r");
rc=EOF+1;
for(i=0; i<b_len && rc!=EOF; i++){
	rc=fscanf(boundary_fp,"%lf,%lf",&boundary_x[i],&boundary_y[i]);	
}
fclose(boundary_fp);
// read in the step size data
rc=EOF+1;
distance_fp = fopen("stepsizes.dat","r");
for(i=0; i<d_len && rc!=EOF; i++){
	rc=fscanf(distance_fp,"%lf",&distance[i]);	
}
fclose(distance_fp);
// open the file to store the distances, loop over step sizes, 
//   calculate the distance, write all of it to file
richdist_fp = fopen("richdist.dat","w");
polypts_fp = fopen("polypts.dat","w");
for(i=0; i<d_len; i++){
	richdist[i]=richardsonDistance(boundary_x,boundary_y,b_len,distance[i],polygon_pts,&polyarea);
	fprintf(richdist_fp,"%.15e, %.15e, %.15e\n",distance[i],richdist[i],polyarea);
	if(i==d_len-1){
		for(j=0; j<2*ceil(richdist[i]/distance[i])+2; j+=2){
			fprintf(polypts_fp,"%.15e, %.15e\n",polygon_pts[j],polygon_pts[j+1]);
		}
	}
}
fclose(richdist_fp);
fclose(polypts_fp);


}
// custom functions
double richardsonDistance(double* X, double* Y, int N, double dist, double* polygon_pts, double* polyArea){
	/* double richardsonDistance(double* X, double* Y, int N, double dist)
	 * ******************************************************************* *
	 * calculates the distance traversed along a path given by (X,Y)
	 * that is measured by lengths of the "ruler" of length "dist."
	 * Input: (X,Y) - arrays representing the tuple (xi,yi) for each index i
	 *        N - length of X and Y
	 *        dist - length of the ruler we are measuring with
	 * Output: d - number of ruler lengths we travelled along the path.
	 *             can be non-integer, since we add the partial ruler 
	 *             length on the last step.
	 * */
	int i=0;
	int j=0;
	int i0=0;
	double d=0.0;
	double dtemp,eps;
	double x,y,m,b,xstar,ystar;
        double polyarea=0.0;

	eps=0.0;
	xstar=X[i0];
	ystar=Y[i0];
	polygon_pts[0]=xstar;
	polygon_pts[1]=ystar;
	j+=2;
	while(i0<N) { //while we haven't reached the end...
		dtemp = 0; //initilize the distance calc variable
		while(dtemp<dist) { // while we haven't reached the "ruler" length
			i++; // increment i 
			x=xstar-X[i]; // calculate the differences in x and y
			y=ystar-Y[i];
			dtemp=sqrt(x*x+y*y); // (2-norm)calculate the distance pythagoras style
			//dtemp=fabs(x)+fabs(y); //(1-norm)
			if(i==N-1){
				d+=dtemp/dist;
				xstar=X[N-1];
				ystar=Y[N-1];
				polygon_pts[j]=xstar;
				polygon_pts[j+1]=ystar;
				polyarea += polygon_pts[j-2]*polygon_pts[j+1]-polygon_pts[j]*polygon_pts[j-1];
				break;
			} 
		}
		if(i==N-1) break;
		//invoke linear interpolation
		linInterp(X[i-1],X[i],Y[i-1],Y[i],&xstar,&ystar,dist);
		if(isnan(xstar)==1||isnan(ystar)==1) { 	
			printf("xstar|ystar = nan, exiting! xstar=%lf,ystar=%lf,dist=%lf,i=%d\n",xstar,ystar,dist,i); 
			break;
		}
		polygon_pts[j]=xstar;
		polygon_pts[j+1]=ystar;
		if(j>=2){
			polyarea += polygon_pts[j-2]*polygon_pts[j+1]-polygon_pts[j]*polygon_pts[j-1];
		}
		// we have gone *exactly* one distance
		d++;// add 1 to the number of ruler lengths travelled
		i0=i; // set the seek index to the next start point
		j+=2;
	}
	(*polyArea)=0.5*fabs(polyarea);
	return d*dist; // return the number of lengths traversed times size of length
}

void linInterp(double x1, double x2,double y1,double y2,double *xstar,double* ystar,double dist){
		int solveInY=0;
		double A,B,C;
	    double m,b,mm,bb,X1,X2,xtemp,ytemp,dx,dy;
	    double disc;
	    dy=y2-y1;
	    dx=x2-x1;
	    if(fabs(dx)<1e-10){// if we are dividing by zero, solve in y
			solveInY=1;
			dx=dy;
			dy=0;
			xtemp=(*xstar);
			(*xstar)=(*ystar);
			(*ystar)=xtemp;
			xtemp=x1;
			x1=y1;
			y1=xtemp;
			xtemp=x2;
			x2=y2;
			y2=xtemp;		
		}
		//linearly interpolate between points to find next place for ruler
		m=(dy)/(dx); //slope
		b=y1-x1*m; //y-intercept
		mm=m*m;//m^2
		bb=b*b;//b^2
		A=1+mm;
		B=2*(m*b-(*xstar)-m*(*ystar));
		C=(*xstar)*(*xstar)+bb-2*b*(*ystar)+(*ystar)*(*ystar)-dist*dist;
		//solve quadratic equation for pythagoras
		disc=B*B-4*A*C;
		disc=(disc<0)?0:sqrt(disc);
		X1=0.5*(-B-disc)/A;
		X2=0.5*(-B+disc)/A;
		//the closest point to this point is the proper point to pick
		(*xstar)=(fabs(x2-X1)<fabs(x2-X2))?X1:X2;
		// find y on the interpolating line
		(*ystar)=m*(*xstar)+b;
		if(solveInY==1){
			xtemp=(*xstar);
			(*xstar)=(*ystar);
			(*ystar)=xtemp;
		}
}

