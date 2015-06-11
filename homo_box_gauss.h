#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//prend en entrÈe img l'image, l'image ou la mettre, z le zoom et n coeff et w,h dimensions

//on fait l'homographie separÈment sur les lignes et les colonnes grace ‡ l'intÈgration
//par convention iimg designe l'image intÈgrale

//fait la somme dans un nouveau vecteur

//En calculant l'écart type de la gaussienne déjà la si on veut un ecart de type de 0.8, on a

// D=2*sqrt(0.64*d^2-0.49);

static int good_modulus_bis(int nn, int p)
{
	if (!p) return 0;
	if (p < 1) return good_modulus(nn, -p);

	unsigned int r;
	if (nn >= 0)
		r = nn % p;
	else {
		unsigned int n = nn;
		r = p - (-n) % p;
		if (r == p)
			r = 0;
	}
	return r;
}


float absf(float a){if(a>0){return a;}else{return -a;}}


float eval_gauss(float *imgw,float j,float D ,int wh){
float sum=0.0; 
float norm=0.0;
int index;
float mean=0.0;
float expo;
if (floor(D)>wh) {
	for (int i=0;i<wh;i++){mean=mean+imgw[0];}
	return mean/((float)(wh)) ;						

	
	}
//if (floor(D)>wh) {return 0;}
for (int i=0;i<D;i++){
	index=floor(j+i-D/2) ;
	//expo=exp(-pow(i-D/2,2)/(2*pow(D/2.0,2)));
	expo=exp(-pow(index-j,2)/(2*pow(D/2.0,2)));
	norm=norm+expo;
	if (floor(j+i-D/2)>0 && floor(j+i-D/2)<wh){
	sum=sum+expo*(imgw[index]);}
}

return  sum/norm   ;}


float gauss_hom_one(float *imgw,float j,float d,int wh){
float a,b,c,dd;
        //j=good_modulus_bis(j,wh);
	float d_aux = 0.64*pow(d,2)-0.49;  //cf. formalisation papier, un peu flou...
	if(d_aux <=0){d_aux=0;}
	float D = 2*sqrt(d_aux);
        //D=1;
	int id=floor(D);
	int ij=floor(j);
	float x = D-id;
	float y = j-ij;
//a l'interieur d'un seul pixel, interpolation bilinaire entre les deux pixels voisin 
	//(il existe mieux que l'interpolation lineaire...))
	/*if(id<=1){a=0 ;
                  b=0 ;
	        if (ij>0 && ij<wh){a=imgw[ij];}
                if (ij+1>0 && ij+1<wh){b=imgw[ij+1];}
		
                  
		return a*(1-y)+b*y;
	}*/
      if(id<=1){a=eval_gauss(imgw,j,1,wh);
		
                  
		return a;
	}
//on est oblige de garantir id petit sinon on sort de iimg3
	//if(id>=wh-1){return 0;}
    
    //ici interpolation bilineaire avec j et d
	a=eval_gauss(imgw,j,D,wh);
	//b=eval_gauss(imgw,j+1,D,wh);
	//c=eval_gauss(imgw,j,D+1,wh);
	//dd=eval_gauss(imgw,j+1,D+1,wh);
	// changer ij et id en j et D pour plus tard...

//(1-x)*(1-y)*a + (1-x)*y*b + x*(1-y)*c + x*y*dd

return a;}







//apply_homo est concu pour H une homographie tel que H[1]=H[4]=H[7]=0
//elle separe la transformation selon les lignes et les colonnes

int apply_homo(float *img,float *img_f,int w,int h,int w_f,int h_f,int mu,int nu,int mu_f,int nu_f,double H[9]){
/*
  * @param
  *     img, img_f : les images d'entrÈe et de sortie
  *     w,h, w_f,h_f : les dimensions des images
  *     mu,nu, mu_f,nu_f : les coordonnÈes du pixel en haut a gauche des images final et initiale
  *     H : homographie telle que b=c=s=0
  * Un pixel ayant une epaisseur de 1, on considere que son antÈcÈdent est d'epaisseur d
  * (d la valeur absolue de la dÈrivÈe de l'homographie en ce point)
  * Dans le code, x et y reprÈsentent les coordonnÈes reelles, float, avec decentrage
  * alors que i et j reprÈsentent les indexes dans le tableau, int, centrÈs en haut ‡ gauche
  * On pourrait Èviter certains dÈcentrage (-mu, -nu) qui seront compensÈs dans linear_int,
  * mais cela permet d'Ítre cohÈrent dans les notations
  */
	int i,j,l;
	
    //w_aux,h_aux, mu_aux,nu_aux pour l'image intermÈdiaire img_aux
    
    int w_aux = w_f; //la 2nde Ètape laisse inchangÈe x, donc w_f=w_aux
    int h_aux = h; //la 1ere Ètape laisse inchangÈe y, donc c'est noir en dehors de cette Èpaisseur
    int mu_aux = mu_f;
    int nu_aux = nu;
	float *img_aux = malloc(w_aux*h_aux*sizeof(float));
	float *img_aux2 = malloc(w_f*h_f*sizeof(float));

    float flmu = (float) mu, flnu_aux = (float) nu_aux;

	for(l=0;l<3;l++){
		float *imgw=malloc(w*sizeof(float)); //une colonne vide

		//operations colonnes par colonnes :
		for(j=0;j<h_aux;j++){
			for(i=0;i<w;i++){imgw[i]=img[(i+w*j)*3+l];}                      	//on extrait le bonne colonne
			
			//build_triple(iimgw,iimgw3,w); 		//on construit la triple integrale, dans une image plus grande

			for(i=0;i<w_aux;i++){
				float x = (float) (i+mu_aux);
				float d = absf((H[0]*H[8]-H[6]*H[2])/pow(H[6]*x+H[8],2)); //derivee selon x
				x = (H[0]*x+H[2])/(H[6]*x+H[8]) - flmu;		//on applique l'homographie
				img_aux[i+j*w_aux] = gauss_hom_one(imgw,x,d,w);		//on realise la convolution
			}
		}
		free(imgw) ; 
		printf("bonjour etape 1\n") ; 

                 //for(i=0;i<w_aux*h_aux;i++){imgw[i]=img[(i+w*j)*3+l]} 

                
		float *imgh=malloc(h_aux*sizeof(float)); //une ligne vide
		
		//opÈrations lignes par lignes, similaire a la precedente :
		for(i=0;i<w_f;i++){
			for(j=0;j<h_aux;j++){imgh[j]=img_aux[i+w_aux*j+l];}  //on extrait le bonne ligne
			
			
			float x =(float) (i+mu_f);
			float d = absf(H[4]/(H[6]*x+H[8])); //dÈrivÈe selon y
			for(j=0;j<h_f;j++){
				float y = (float) (j+nu_f);
				y = (H[4]*y+H[5])/(H[6]*x+H[8]) - flnu_aux;
				img_aux2[i+j*w_f] = gauss_hom_one(imgh,y,d,h_aux);
			}
		}
		for(i=0;i<w_f;i++){
                       for(j=0;j<h_f;j++)   {img_f[3*(i+w_f*j)+l]=img_aux2[i+w_f*j];}}

	free(imgh);
		}
                
	return 0;
}
