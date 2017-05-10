/**Zoetgnandé Yannick Wend Kuni*/



/**Inclusion des bibliothèques dont on aura bésoin*/
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <math.h>
#include <time.h>

/**Definition des constantes*/
#define maxD 64
#define maxT 4096
#define pixelMaxHauteur 1200/*La hauteur maximale pour une image*/
#define pixelMaxLargeur 1600/*La largeure maximale pour une image*/
#define imageMaxTab 1920000/*La taille maximale du tableau temporaire qui va contenir les bits
                            avant d'être mis dans la matrice de traitement*/
#define parts 8/*La taille pour chaque sous matrice 8*8*/


/**Definition des types de variables personnalisées*/
typedef double matD[maxD][maxD];
typedef double tabD[maxT];
typedef unsigned char imageMatrice[pixelMaxLargeur][pixelMaxHauteur];/*Matrice qui va contenir la matrice de l'image*/
typedef unsigned char imageTableau[imageMaxTab];/*Tableau qui va contenir les bits de la matrice image*/
typedef unsigned char matrice8x8[parts][parts];/*Matrice qui represente les matrices de tailles 8*8*/
typedef double matrice8x8D[parts][parts];/*Matrice DCT et IDCT (double) qui represente les matrices de tailles 8*8*/
typedef double tableau64[64];

/*Une liste chainée qui va contenir les matrices 8*8*/
typedef struct noeud
{
    struct noeud *suiv;
    matrice8x8 M;
};

/*Une liste chainee qui va contenir les matrice doubles 8*8*/
typedef struct node
{
    struct node *suiv;
    matrice8x8D M;
};

/*Une liste qui va contenir les tableaux issus du RLE tout juste après le parcours Zigzag de la liste chainee*/
typedef struct tabM8X8
{
    struct tabM8X8 *suiv;
    tableau64 T;
};


typedef noeud* liste;

typedef node* liste_double;

typedef tabM8X8* parcours;

double cx(int x)
{
    if(x>0)
        return 1;
    else
        return (1.0/sqrt(2));
}
void RLE_compresse(tabD T, int n, tabD TC, int &nc)
{
    int i,j,k,l,nbr,nb;
    j=0;
    nbr=1;
    nb=0;
    for(i=1;i<=n;i++)
    {
        if(T[i-1]==T[i])
            {
                nbr++;
                if(i==n)
                {
                    TC[j++]=nbr;
                    TC[j++]=T[i-1];
                    nc=j;
                }
            }
        else
            {
                TC[j++]=nbr;
                TC[j++]=T[i-1];
                nc=j;
                nbr=1;
            }
    }
}
void RLE_decompresse(tabD TC, int nc, tabD TI, int &n)
{
    int i,j,k;
    n=0;
    i=0;
    while(i<nc-1)
    {
        for(j=0;j<TC[i];j++)
        {
            TI[n]=TC[i+1];
            n++;
        }
        i=i+2;
    }
}
void DCT_matrice(matrice8x8 M, matrice8x8D D,int n)
{
    int i,j,x,y;
    double res=0;
    const double pi=22.0/7.0;
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            for(x=0;x<n;x++)
            {
                for(y=0;y<n;y++)
                {
                    res+=M[x][y]*cos(((2*x+1)*i*pi)/(2.0*n))*cos(((2*y+1)*j*pi)/(2.0*n));
                }
            }
            D[i][j]=(2.0/n)*cx(i)*cx(j)*res;
            //break;
            res=0;
        }
}
void IDCT_matrice(matrice8x8D D, matrice8x8 IDCT, int n)
{
    int i,j,x,y;
    double res;
    const double pi=22.0/7.0;
    for(x=0;x<n;x++)
        for(y=0;y<n;y++)
        {
            res=0;
            for(i=0;i<n;i++)
            {
                for(j=0;j<n;j++)
                {
                    res+=cx(i)*cx(j)*D[i][j]*cos(((2*x+1)*i*pi)/(2.0*n))*cos(((2*y+1)*j*pi)/(2.0*n));
                }
            }
            IDCT[x][y]=(2.0/n)*res;
        }
}
void remplirMatriceAleatoire(matrice8x8 M,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            M[i][j]=((rand() / (double)RAND_MAX * (9)));
        }
    }
}
void afficherMatrice(matrice8x8 M,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%2.0f\t",M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
void afficherMatriceDouble(matrice8x8D M,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%2.0f\t",M[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
int matrice_egale(matrice8x8 M1, matrice8x8 M2,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
           if(fabs(M1[i][j]-M2[i][j])<0.001)
                return 0;
        }
    }
    return 1;
}
int difference(matrice8x8 M1,matrice8x8 M2,matrice8x8 M3,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            M3[i][j]=M1[i][j]-M2[i][j];
        }
    }
}
void CalculMatriceQuant(matrice8x8D Q,int n)
{
    int k,i,j;
    printf("donner un facteur de qualité entre 1 et 25\n");
    scanf("%d",&k);
    while((k<0) || (k>25))
    {
        printf("donner un facteur de qualité entre 1 et 25\n");
        scanf("%d",&k);
    }
    printf("\nk=%d\n\n",k);
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            {
                Q[i][j]=1+k*(1+i+j);
            }
    }
}
void Quantification(matrice8x8D D,matrice8x8D DQ,matrice8x8D Q,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            {
                DQ[i][j]=D[i][j]/Q[i][j];
            }
    }
}
void Dequantification(matrice8x8 D,matrice8x8 DQ,matrice8x8 Q,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
            {
                DQ[i][j]=D[i][j]*Q[i][j];
            }
    }
}
void Zigzag(matrice8x8D D,int n,tableau64 Z)
{
    int i,j,k,d,x,y;
    i=0;
    j=0;
    k=0;
        for(d=0;d<n;d++)
        {
            if(d%2==0)
            {
                for(i=d,j=0;i>=0,j<=d;i--,j++)
                {
                    Z[k]=D[i][j];
                    k++;
                }
            }
            if(d%2!=0)
            {
                for(j=d,i=0;i<=d,j>=0;i++,j--)
                {
                    Z[k]=D[i][j];
                    k++;
                }
            }
        }
        while(d<(2*n-1))
        {
            if(d%2==0)
            {
                for(i=n-1,j=d-n+1;i>=d-n+1,j<=n-1;i--,j++)
                {
                    Z[k]=D[i][j];
                    k++;
                }
            }
            if(d%2!=0)
            {
                for(j=n-1,i=d-n+1;i<=n-1,j>=d-n+1;i++,j--)
                {
                    Z[k]=D[i][j];
                    k++;
                }
            }
            d++;
        }
}
void RetourZigzag(matD D,int n,tabD Z)
{
    int i,j,k,d,x,y;
    i=0;
    j=0;
    k=0;
        for(d=0;d<n;d++)
        {
            if(d%2==0)
            {
                for(i=d,j=0;i>=0,j<=d;i--,j++)
                {
                    D[i][j]=Z[k];
                    k++;
                }
            }
            if(d%2!=0)
            {
                for(j=d,i=0;i<=d,j>=0;i++,j--)
                {
                    D[i][j]=Z[k];
                    k++;
                }
            }
        }
        while(d<(2*n-1))
        {
            if(d%2==0)
            {
                for(i=n-1,j=d-n+1;i>=d-n+1,j<=n-1;i--,j++)
                {
                    D[i][j]=Z[k];
                    k++;
                }
            }
            if(d%2!=0)
            {
                for(j=n-1,i=d-n+1;i<=n-1,j>=d-n+1;i++,j--)
                {
                    D[i][j]=Z[k];
                    k++;
                }
            }
            d++;
        }
}
void afficherTableau(tabD T,int n)
{
    int i,j;
    for(i=0;i<n;i++)
    {
            printf("Tab[%d]=%2.0f\n",i,T[i]);
    }
}
void InsererMatrice(liste &l,matrice8x8 M8)
{
    int i,j;
    /* On crée un nouvel élément */
    noeud* q;
    q =(noeud*) malloc(sizeof(noeud));
    /* On assigne la valeur au nouvel élément */
    if(l==NULL)
    {
        for(i=0;i<8;i++)
        {
            for(j=0;j<8;j++)
            {
                q->M[i][j]=M8[i][j];
            }
        }
        q->suiv=NULL;
        l=q;
    }
    else
    {
        liste p=l;
        while(p->suiv!=NULL)
        {
            p=p->suiv;
        }
        for(i=0;i<8;i++)
        {
            for(j=0;j<8;j++)
            {
                q->M[i][j]=M8[i][j];
            }
        }
        q->suiv=NULL;
        p->suiv=q;
    }
}
void InsererMatriceDouble(liste_double &l,matrice8x8D M8)
{
    int i,j;
    /* On crée un nouvel élément */
    node* q;
    q =(node*) malloc(sizeof(node));
    /* On assigne la valeur au nouvel élément */
    if(l==NULL)
    {
        for(i=0;i<8;i++)
        {
            for(j=0;j<8;j++)
            {
                q->M[i][j]=M8[i][j];
            }
        }
        q->suiv=NULL;
        l=q;
    }
    else
    {
        liste_double p=l;
        while(p->suiv!=NULL)
        {
            p=p->suiv;
        }
        for(i=0;i<8;i++)
        {
            for(j=0;j<8;j++)
            {
                q->M[i][j]=M8[i][j];
            }
        }
        q->suiv=NULL;
        p->suiv=q;
    }
}
void diviserMatrice(imageMatrice M,int n,int m,liste &l)
{
    int i,j,k,p,s,r,a,b;
    matrice8x8 M8;
    s=n/8;
    r=m/8;
    k=0;
    p=0;
    while(k<s)
    {
        while(p<r)
        {
            for(i=0;i<8;i++)
            {
                a=i+8*(k);
                for(j=0;j<8;j++)
                {
                    b=j+8*(p);
                    M8[i][j]=M[a][b];
                    //printf("\na=%d\tb=%d\n",a,b);
                    //printf("%c\t",M8[i][j]);
                }
                //printf("\n");
            }
            //printf("\n");
            InsererMatrice(l,M8);
            p++;
        }
        //printf("\n\n");
        k++;
        p=0;
    }
}
void afficherListe(liste l)
{
    int i,j,k=0;
    liste* p;
    //FILE * f;
    //f = fopen("resultat.txt", "w");
    if(l!=NULL)
    {
        while(l!=NULL)
        {
            printf("Affichage de la %d eme matrice:\n",k);
            for(i=0;i<8;i++)
            {
                for(j=0;j<8;j++)
                {
                    printf("%c\t",l->M[i][j]);
                }
                printf("\n");
            }
            l=l->suiv;
            k++;
            //putc('\n',f);
            //putc('\n',f);
            printf("\n\n\n");
        }
    }
    else
    {
        printf("La liste est vide\n");
    }
    //fclose(f);

}
void afficherListeDouble(liste_double l)
{
    int i,j,k=0;
    liste_double* p;
    //FILE * f;
    //f = fopen("resultat.txt", "w");
    if(l!=NULL)
    {
        while(l!=NULL)
        {
            printf("Affichage de la %d eme matrice:\n",k);
            for(i=0;i<8;i++)
            {
                for(j=0;j<8;j++)
                {
                    printf("%2.0f\t",l->M[i][j]);
                }
                printf("\n");
            }
            l=l->suiv;
            k++;
            //putc('\n',f);
            //putc('\n',f);
            printf("\n\n\n");
        }
    }
    else
    {
        printf("La liste est vide\n");
    }
    //fclose(f);

}
void DCT_liste_matrice(liste l,liste_double &ld)
{
    int i,j,n;
    n=8;
    matrice8x8D D8;
    if(l!=NULL)
    {
        while(l!=NULL)
        {
            DCT_matrice(l->M,D8,n);
            //InsererMatrice(ld,D8);
            InsererMatriceDouble(ld,D8);
            l=l->suiv;
            //printf("\n\n\n");

        }
    }
    else
    {
        printf("\nLa liste est vide\n");
    }
}
void QuantificationListe(matrice8x8D Q8,liste_double l,liste_double &ld)
{
    int i,j,n=8;
    matrice8x8D DQ8;
    if(l!=NULL)
    {
        while(l!=NULL)
        {
            Quantification(l->M,DQ8,Q8,n);
            InsererMatriceDouble(ld,DQ8);
            l=l->suiv;
        }
    }
    else
    {
        printf("\nLa liste est vide\n");
    }
}
void InsererTableauDouble(parcours &p,tableau64 tab)
{
    int i,j;
    /* On crée un nouvel élément */
    tabM8X8* q;
    q =(tabM8X8*) malloc(sizeof(tabM8X8));
    /* On assigne la valeur au nouvel élément */
    if(p==NULL)
    {
        for(i=0;i<64;i++)
        {
                q->T[i]=tab[i];
                //printf("\n%f\n",q->T[i]);
        }
        q->suiv=NULL;
        p=q;
    }
    else
    {
        parcours l=p;
        while(l->suiv!=NULL)
        {
            l=l->suiv;
        }
        for(i=0;i<64;i++)
        {
                q->T[i]=tab[i];
        }
        q->suiv=NULL;
        l->suiv=q;
    }
}
void ParcoursZigzagTout(liste_double l,parcours &p)
{
    int i,j;
    tableau64 tab;
    if(l!=NULL)
    {
        while(l!=NULL)
        {
            Zigzag(l->M,8,tab);
            InsererTableauDouble(p,tab);
            l=l->suiv;
        }
    }
    else
    {
        printf("\nLa liste est vide\n");
    }
}
void afficherTableauListe(parcours p)
{
    int i,k=1;
    if(p==NULL)
    {
        printf("\nListe vide\n");
    }
    else
    {
        while(p!=NULL)
        {
            printf("\nAffichage du %d eme Tableau:\n",k);
            for(i=0;i<64;i++)
            {
                printf("%2.0f\t",p->T[i]);
            }
            p=p->suiv;
        }
    }
}
main()
{
    int n,nt,nc,nd,i,j,k,longueur,largeur,m;
    matD D,M,IDCT,DQ,RZ,DQuant,MI;
    tabD Z,ZC,DZ;
    matrice8x8D Q;
    n=8;
    imageMatrice Mat;
    int entete[10];
    int r=0,s=0;
    char src[]="lena_bw_p1.pbm";
    FILE * fsrc;
    FILE * f;
    fsrc = fopen(src, "r");
    f = fopen("hello.txt", "w");
    if (fsrc == NULL)
    {
        perror(src);
    }
    else
    {
                int c;
                i=0;
                fscanf(fsrc, "%s %d %d", &entete[0], &entete[1], &entete[2]);
                longueur=entete[1];
                largeur=entete[2];
                m=longueur*largeur;
                imageTableau Tab;
                while ((c = getc(fsrc)) != EOF)
                {
                    if((c=='0') || (c=='1'))
                    {
                        Tab[i]=c;
                        i++;
                    }
                    r++;
                }
                if(i==m)
                {
                    i=0;
                    k=0;
                    for(i=0;i<longueur;i++)
                    {
                        for(j=0;j<largeur;j++)
                        {
                            Mat[i][j]=Tab[k];
                            //putc(Mat[i][j],f);
                            //putc('\n',f);
                            k++;
                        }
                    }
                }
    }
    liste l=NULL;
    liste_double liste_DCT=NULL,liste_quant=NULL;
    parcours p=NULL;
    diviserMatrice(Mat,512,512,l);
    DCT_liste_matrice(l,liste_DCT);
    afficherListeDouble(liste_DCT);
    CalculMatriceQuant(Q,n);
    printf("Affichage de la matrice de quantification\n");
    afficherMatriceDouble(Q,n);
    QuantificationListe(Q,liste_DCT,liste_quant);
    afficherListeDouble(liste_quant);
    ParcoursZigzagTout(liste_quant,p);
    afficherTableauListe(p);
    fclose(fsrc);
}
