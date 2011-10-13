#include <R.h>
#include <Rmath.h>

void isbfC(double *residu, double *seuille, double *beta, double *s, double *impmin, int *p, int *K)
{

   // Déclaration des variables
   int i,i1,j,colonne,col1,ibest,jbest,colbest,ib1,vmin,vmax;
   double calcul, imp=(*impmin)+1, signe;

   // Remplissage des tableaux "résidus" (contenant les sommes des blocs) et "seuille"
   for (j=0;j<(*p);j++)
   {
      calcul = fabs(residu[j])-(*s);
      if (calcul>0) seuille[j] = calcul*calcul;
   }
   for (i=1;i<(*K);i++)
   {
      i1 = i+1;
      col1 = (i-1)*(*p);
      colonne = col1+(*p);
      for (j=0;j<(*p)-i;j++)
      {
         residu[colonne+j] = residu[col1+j] + residu[j+i];
         calcul = fabs(residu[colonne+j])/(i1)-(*s)/sqrt(i1);
         if (calcul>0) seuille[colonne+j] = calcul*calcul*(i1);
      }
   }
   
   // Recherche du bloc impliquant la meilleure amélioration
   while (imp>(*impmin))
   {
      imp = 0;
      ibest = 0;
      colbest = 0;
      for (j=0;j<(*p);j++) if (seuille[j]>imp)
      {
         jbest = j;
         imp = seuille[j];
      }
      for (i=1;i<(*K);i++)
      {
         colonne = i*(*p);
         for (j=0;j<(*p)-i;j++) if (seuille[colonne+j]>imp)
         {
            ibest = i;
            jbest = j;
            colbest = colonne;
            imp = seuille[colonne+j];
         }
      }

      // Mise à jour, en conséquence, de "beta", et de "résidu" et "seuille"
      if (imp>(*impmin))
      {
         ib1 = ibest +1;
         calcul = residu[colbest+jbest]/ib1;
         if (calcul>0) signe = 1; else signe = -1;
         calcul = signe*(fabs(calcul)-(*s)/sqrt(ib1));
         for (j=jbest;j<jbest+ib1;j++)
         {
            beta[j] = beta[j] + calcul;
            residu[j] = residu[j] - calcul;
            signe = fabs(residu[j])-(*s);
            if (signe>0) seuille[j] = signe*signe; else seuille[j] = 0;
         }
         for (i=1;i<(*K);i++)
         {
            i1 = i+1;
            col1 = (i-1)*(*p);
            colonne = col1+(*p);
            vmin = jbest-i;
            if (vmin<0) vmin = 0;
            vmax = jbest + ib1;
            if (vmax>(*p)-i) vmax = (*p)-i;
            for (j=vmin;j<vmax;j++)
            {
               residu[colonne+j] = residu[col1+j] + residu[j+i];
               calcul = fabs(residu[colonne+j])/(i1)-(*s)/sqrt(i1);
               if (calcul>0) seuille[colonne+j] = calcul*calcul*(i1); else seuille[colonne+j] = 0;
            }
         }
      }
   }

   // Fin de la fonction
   return;
}

void isbfRegC(double *X, double *Y, double *COV, double *beta, double *s, double *impmin, double *fgroups, int *p, int *n, int *K)
{

   // Déclaration des variables
   int h, i, j, l, cpt, cpt2, ibest, hbest;
   double imp=(*impmin)+1, S, S2, coef, cov, calcul, coefbest;

   // Calcul préliminaire des 'Covariances' (renormalisations)
   cpt = 0;
   for (i=1;i<(*K)+1;i++)
   {
      cpt2 = (*p)-i+1;
      for (h=0;h<cpt2;h++)
      {
         S = 0;
         for (l=0;l<(*n);l++)
         {
            S2 = 0;
            for (j=0;j<i;j++) S2 += X[(*n)*(h+j)+l];
            S += S2*S2;
         }
         COV[cpt+h] = S;
      }
      cpt += cpt2;
   }

   // Recherche du bloc impliquant la meilleure amélioration
   while (imp>(*impmin))
   {
      imp = 0;
      cpt = 0;
      for (i=1;i<(*K)+1;i++)
      {
         cpt2 = (*p)-i+1;
         for (h=0;h<cpt2;h++)
         {
            cov = COV[cpt+h];
            S = 0;
            for (l=0;l<(*n);l++)
            {
               S2 = 0;
               for (j=0;j<i;j++) S2 += X[(*n)*(h+j)+l];
               S += S2*Y[l];
            }
            coef = S/cov;
            calcul = (1+(*fgroups)/sqrt(i))*(*s)/sqrt(cov);
            if (coef>calcul) coef -= calcul; else if (coef<(-calcul)) coef += calcul; else coef = 0;
            calcul = coef*coef*cov;
            if (calcul>imp)
            {
               ibest = i;
               hbest = h;
               imp = calcul;
               coefbest = coef;
            }
         }
         cpt += cpt2;
      }

      // Mise à jour, en conséquence, de "beta", et de "Y"
      if (imp>(*impmin))
      {
         for (j=0;j<ibest;j++) beta[hbest+j] += coefbest;
         for (l=0;l<(*n);l++)
         {
            S2 = 0;
            for (j=0;j<ibest;j++) S2 += X[(*n)*(hbest+j)+l]*coefbest;
            Y[l] -= S2;
         }
      }
   }

   // Fin de la fonction
   return;
}

