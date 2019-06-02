/*
 * fft.cpp
 *
 *  Created on: May 14, 2019
 *      Author: vamshi.kurva
 */

#include "opencv/highgui.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv/cv.h"
#include <iostream>
#include <complex>
#include <cmath>
using namespace cv;
using namespace std;
class AU_FFT
{
	typedef complex<double> comp;
    public:
	Mat img,mag,phase,inv;
    int REF,SIZE;
    comp** out;
    class rever
    {
     public:
     int a,b;
    };
	 AU_FFT(Mat img1,int a)       //constructor
	 {
		  int i,j;
		  SIZE=a;
		  img.create(SIZE,SIZE,CV_8UC1);
		  int m,n,M,N;
		  M=img1.rows;
		  N=img1.cols;
		  //cout<<" "<<M<<" "<<N;
		  m=ceil(log2(M));
		  n=ceil(log2(N));
		  M=pow(2,m);
		  N=pow(2,n);
		  if(M>N)
		    N=M;
		  else
			M=N;
		  for(m=0;m<M;m++)
		  {
			  for(n=0;n<M;n++)
			  {
                // extending the image for continuity
				  /*
				 i=m/img1.rows;
				 j=n/img1.cols;
				 if((3*i)%2 ==0 && (3*j)%2==0)
					 img.at<unsigned char>(m,n)=img1.at<unsigned char>(m-i*img1.rows,n-j*img1.cols);
				 else if((3*i)%2==0 && (3*j)%2!=0)
					 img.at<unsigned char>(m,n)=img1.at<unsigned char>(m-i*img1.rows,img1.cols*(1+j)-n-1);
				 else if((3*i)%2!=0 && (3*j)%2==0)
				 	 img.at<unsigned char>(m,n)=img1.at<unsigned char>(img1.rows*(1+i)-m-1,n-j*img1.cols);
				 else
				 	 img.at<unsigned char>(m,n)=img1.at<unsigned char>(img1.rows*(1+i)-m-1,img1.cols*(1+j)-n-1);
				 	 */
			   // padding zeros
				   if(m<img1.rows && n<img1.cols)
						img.at<unsigned char>(m,n)=img1.at<unsigned char>(m,n);
					 else
						img.at<unsigned char>(m,n)=0;
			  }
		  }
		 out=new comp*[SIZE];
         for(m=0;m<SIZE;m++)
        	 out[m]=new comp[SIZE];
		 mag=img.clone();
		 phase=img.clone();
		 inv=img.clone();
	 }
     void dft1(comp** ,comp** ,int ,int ,int ,int ,char);
     void fft(int );
     void shift();
     void ifft(int );
     ~AU_FFT()         //destructor
     {
    	 int m;
    	 for(m=0;m<SIZE;m++)
    	 {
    		 delete[] out[m];
    	 }
    	 delete[] out;
     }
};
void AU_FFT::dft1(comp** ptr1,comp** ptr2,int pos1,int pos2,int fact,int N,char e)
{
 int k=1,l1,a,b,t1,t2,p=1,m,n,M1,N1,M;
 M1=N1=M=N;
 comp static x[4],temp1,temp2,j,**p1;
 comp **p3,**p2;
 p3=new comp*[SIZE];
 p2=new comp*[SIZE];
 p1=new comp*[SIZE];
 for(m=0;m<SIZE;m++)
 {
	 p3[m]=new comp[SIZE];
	 p2[m]=new comp[SIZE];
	 p1[m]=new comp[SIZE];
 }
 for(m=0;m<SIZE;m++)
 {
	 for(n=0;n<SIZE;n++)
	 {
		 p3[m][n]=ptr1[m][n];
		 p2[m][n]=ptr2[m][n];
	 }
 }
 j=-1;
 j=sqrt(j);
 double pi,w1,w2;
 pi=2*asin(1);
 for(m=0;m<N;m++)
 {
   for(n=0;n<N;n++)
    {
      p1[m][n]=p3[m+pos1][n+pos2];
    }
  }
 if(e=='D')
	 j=j;
 if(e=='I')
	 j=-j;
 //direct method for 2D-FFT
 do
 {
   l1=m=t1=0;
   for(a=0;a<M/2;a++)
   {
   	if(m==l1+M/pow(2,k))
   	{
   	    if(a==M/2)
   	         break;
            m=m+M/pow(2,k);
            l1=m;
            a--;
            continue;
    }
         n=t2=0;
         int l2=0;
         for(b=0;b<N/2;b++)
          {
            if(n==l2+N/pow(2,p))
   	        {
   	           if(b==N/2)
   	             break;
               n=n+N/pow(2,p);
               l2=n;
               b--;
               continue;
            }
            w1=(2*pi)/(M1);
            w2=(2*pi)/(N1);
            temp1=exp(-(w1*t1*j));
            temp2=exp(-(w2*t2*j));
            	x[0]=p1[m][n]+p1[m+M1/2][n]+p1[m][n+N1/2]+p1[m+M1/2][n+N1/2];
            	x[1]=(p1[m][n]+p1[m+M1/2][n]-p1[m][n+N1/2]-p1[m+M1/2][n+N1/2])*temp2;
            	x[2]=(p1[m][n]-p1[m+M1/2][n]+p1[m][n+N1/2]-p1[m+M1/2][n+N1/2])*temp1;
            	x[3]=(p1[m][n]-p1[m+M1/2][n]-p1[m][n+N1/2]+p1[m+M1/2][n+N1/2])*temp1*temp2;
           p1[m][n]=x[0];
           p1[m][n+N1/2]=x[1];
           p1[m+M1/2][n]=x[2];
           p1[m+M1/2][n+N1/2]=x[3];
           n++;
           t2++;
           if(t2==(N1/2))
        	   t2=0;
        }
         m++;
         t1++;
         if(t1==(M1/2))
        	 t1=0;
   }
  M1=M1/2;
  N1=N1/2;
  k++;
  p++;
 }while(k<=log2(M));
 // bit reversal order
 // cout<<endl<<"fft is....";
 int *ptr=new int[N];
 rever y[N][N];
 int i,l;
 // 1D-bit reversal
 n=log2(N);
 ptr[0]=p=0;
 ptr[1]=1;
 for(i=1;i<n;i++)
 {
    p++;
    l=pow(2,p);
    for(m=0;m<l;m++)
    {
       ptr[m]=2*ptr[m];
       ptr[m+l]=ptr[m]+1;
    }
 }
 // 2D-bit reversal
 for(i=0;i<N;i++)
 {
   for(m=0;m<N;m++)
   {
    y[i][m].a=ptr[i];
    y[i][m].b=ptr[m];
   }
 }
 for(m=0;m<N;m++)
 {
   for(n=0;n<N;n++)
    {
    //  p2[m+pos1][n+pos2]=p1[y[m][n].a][y[m][n].b];
      ptr2[m+pos1][n+pos2]=p1[y[m][n].a][y[m][n].b];
    }
 }
 delete ptr;
 for(m=0;m<SIZE;m++)
 {
	 delete[] p3[m];
	 delete[] p2[m];
	 delete[] p1[m];
 }
 delete[] p3;
 delete[] p2;
 delete[] p1;
}

void AU_FFT::fft(int a)
{
	// cout<<" vamshi";
	 comp  **p1,**p2,**p3,im,sum;
	 int fact,M,N,m,n,i,j,i1,j1,i2,j2,k1,k2;
	 p3=new comp*[SIZE];
	 p2=new comp*[SIZE];
	 p1=new comp*[SIZE];
	 for(m=0;m<SIZE;m++)
	 {
		 p3[m]=new comp[SIZE];
		 p2[m]=new comp[SIZE];
		 p1[m]=new comp[SIZE];
	 }
	 double pi,w,mul;
	 char e;
	 REF=a;
	 M=img.rows;
	 N=img.cols;
	 for(i=0;i<img.rows;i++)
	 {
	   for(j=0;j<img.cols;j++)
	    {
	       p1[i][j].real(img.at<unsigned char>(i,j));
	       p1[i][j].imag(0);
	    }
	 }
	 for(i=img.rows;i<SIZE;i++)
	 {
		for(j=img.cols;j<SIZE;j++)
		{
			p1[i][j]=0;
		}
	 }
	  pi=2*asin(1);
	  im=-1;
	  im=sqrt(im);
	  m=ceil(log2(M));
	  M=pow(2,m);
	  n=ceil(log2(N));
	  N=pow(2,n);
	  if(M>N)
	    N=M;
	  else
		M=N;

	  mul=sqrt(M*M);
	  mul=1/mul;
	  fact=M/REF;
	  w=2*pi/M;
	  e='D';
	  for(i1=0;i1<fact;i1++)
	  {
	    for(j1=0;j1<fact;j1++)
	     {
	        for(m=0;m<M;m++)
	        {
	          for(n=0;n<M;n++)
	           {
	             p2[m][n]=p1[m][n];
	           }
	        }
	        for(m=0;m<M;m++)
	        {
	          for(n=0;n<M;n++)
	           {
	              p2[m][n]=p2[m][n]*exp(-(w*m*i1*im))*exp(-(w*n*j1*im));
	           }
	        }
	        for(i2=0;i2<fact;i2++)
	        {
	          for(j2=0;j2<fact;j2++)
	           {
	        	  dft1(p2,p3,i2*REF,j2*REF,fact,REF,e);
	           }
	        }
	       for(k1=0;k1<REF;k1++)
	       {
	         for(k2=0;k2<REF;k2++)
	          {
	            sum=0;
	            for(m=0;m<fact;m++)
	            {
	              for(n=0;n<fact;n++)
	               {
	                 sum=sum+p3[k1+m*REF][k2+n*REF];
	               }
	            }
	            out[fact*k1+i1][fact*k2+j1]=sum;
	          }
	       }
	    }
	  }
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		 {
			 out[i][j]*=mul;
			 w=(abs(out[i][j])*255*100)/abs(out[0][0]);
			 phase.at<unsigned char>(i,j)=arg(out[i][j]);
			// cout<<" "<<int(phase.at<unsigned char>(i,j));
		   	 mag.at<unsigned char>(i,j)=(abs(out[i][j])); //w;
		 }
	}
	//cout<<" "<<"everything is fine";
	 for(m=0;m<SIZE;m++)
	 {
		 delete[] p3[m];
		 delete[] p2[m];
		 delete[] p1[m];
	 }
	 delete[] p3;
	 delete[] p2;
	 delete[] p1;
}
void AU_FFT::shift()
{
	int i,j,M,N;
	M=mag.rows;
	N=mag.cols;
    for(i=0;i<M/2-1;i++)
    {
      for(j=0;j<N/2-1;j++)
      {
    	   //original spectrum
    	  mag.at<unsigned char>(i,j)=abs(out[i+M/2][j+M/2]);
    	  mag.at<unsigned char>(i+M/2,j+M/2)=abs(out[i][j]);
    	  mag.at<unsigned char>(i+M/2,j)=abs(out[i][j+M/2]);
    	  mag.at<unsigned char>(i,j+M/2)=abs(out[i+M/2][j]);

    	  phase.at<unsigned char>(i,j)=arg(out[i+M/2][j+M/2]);
    	  phase.at<unsigned char>(i+M/2,j+M/2)=arg(out[i][j]);
    	  phase.at<unsigned char>(i+M/2,j)=arg(out[i][j+M/2]);
    	  phase.at<unsigned char>(i,j+M/2)=arg(out[i+M/2][j]);
      }
    }
}
void AU_FFT::ifft(int a)
{
	 comp **p2,**p3,**p1,im,sum;
	 int fact,M,N,m,n,i,j,i1,j1,i2,j2,k1,k2;
	 p3=new comp*[SIZE];
	 p2=new comp*[SIZE];
	 p1=new comp*[SIZE];
	 for(m=0;m<SIZE;m++)
	 {
		 p3[m]=new comp[SIZE];
		 p2[m]=new comp[SIZE];
		 p1[m]=new comp[SIZE];
	 }
	 double pi,w,mul;
	 char e;
	 im=-1;
	 REF=a;
	 im=sqrt(im);
	 M=mag.rows;
	 N=mag.cols;
	 //if(M==N)
		// cout<<" ok "<<M;
	  pi=2*asin(1);
	  mul=sqrt(M*M);
	  mul=1/mul;
	  fact=M/REF;
	  w=2*pi/M;
	  e='I';
	  for(i1=0;i1<fact;i1++)
	  {
	    for(j1=0;j1<fact;j1++)
	     {
	        for(m=0;m<M;m++)
	        {
	          for(n=0;n<N;n++)
	           {
	             p2[m][n]=out[m][n];
	           }
	        }
	        for(m=0;m<M;m++)
	        {
	          for(n=0;n<N;n++)
	           {
	              p2[m][n]=p2[m][n]*exp((w*m*i1*im))*exp((w*n*j1*im));
	           }
	        }
	        for(i2=0;i2<fact;i2++)
	        {
	          for(j2=0;j2<fact;j2++)
	           {
	        	  dft1(p2,p3,i2*REF,j2*REF,fact,REF,e);
	           }
	        }
	       for(k1=0;k1<REF;k1++)
	       {
	         for(k2=0;k2<REF;k2++)
	          {
	            sum=0;
	            for(m=0;m<fact;m++)
	            {
	              for(n=0;n<fact;n++)
	               {
	                 sum=sum+p3[k1+m*REF][k2+n*REF];
	               }
	            }
	            p1[fact*k1+i1][fact*k2+j1]=sum;
	          }
	       }
	    }
	  }
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		 {
			 p1[i][j]*=mul;
		   	 inv.at<unsigned char>(i,j)=abs(p1[i][j]);
		 }
	}
	//cout<<" "<<"everything is fine";
	 for(m=0;m<SIZE;m++)
	 {
		 delete[] p3[m];
		 delete[] p2[m];
		 delete[] p1[m];
	 }
	 delete[] p3;
	 delete[] p2;
	 delete[] p1;
}
int main()
{
	  Mat img1 = imread("test61.jpeg", 0 );
	  if (img1.empty()) //1024 check whether the image is loaded or not
	  {
	  	cout << "Error : Image cannot be loaded..!!" << endl;
	  	       //system("pause"); //wait for a key press
	  	return -1;
	  }
	  int m,n,M,N;
	  M=img1.rows;
	  N=img1.cols;
	 // cout<<" "<<M<<" "<<N;
	  m=ceil(log2(M));
	  n=ceil(log2(N));
	  M=pow(2,m);
	  N=pow(2,n);
	  if(M>N)
	    N=M;
	  else
		M=N;

	AU_FFT IMG(img1,M);
	namedWindow("original", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imwrite("padded_image.png", IMG.img);
	imshow("original", IMG.img); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("original");
   // imwrite("pad2.jpeg",IMG.img);

    IMG.fft(M);
	namedWindow("Magnitude spectrum", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imwrite("Mag_spectrum.png", IMG.mag);
	imshow("Magnitude spectrum", IMG.mag); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("Magnitude spectrum");

	namedWindow("phase spectrum", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imshow("phase spectrum", IMG.phase); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("phase spectrum");

    IMG.shift();
	namedWindow("shifted mag spectrum", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imwrite("shifted_mag_spectrum.png", IMG.mag);
	imshow("shifted mag spectrum", IMG.mag); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("shifted mag spectrum");

	namedWindow("shifted phase spectrum", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imshow("shifted phase spectrum", IMG.phase); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("shifted phase spectrum");

    //inverse
    IMG.ifft(M);
	namedWindow("inverse spectrum", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imshow("inverse spectrum", IMG.inv); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
    destroyWindow("inverse spectrum");

     //difference image
	 for(int i=0;i<img1.rows;i++)
	 {
	 	for(int j=0;j<img1.cols;j++)
	 	{
	 		img1.at<unsigned char>(i,j)=img1.at<unsigned char>(i,j)-IMG.inv.at<unsigned char>(i,j);
	 	 //   error=error+pow((fft1.at<unsigned char>(i,j)),2);
	 	}
	 }
	namedWindow("difference image", CV_WINDOW_AUTOSIZE); //create a window with the name "MyWindow"
	imshow("difference image", img1); //display the image which is stored in the 'img' in the "MyWindow" window
	waitKey(0); //wait infinite time for a keypress
	destroyWindow("difference image");
	//end */

  return 0;
}
