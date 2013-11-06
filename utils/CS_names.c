
/* 
   This fragment  should be added to the file niiter.c of chilespice.
   There should be a similar fragment already in the chilespice repository, but
   without the printouts to the file namesChile.txt.  This code should replace
   that older fragment.


   What this fragment does, is add the capability for chilespice to output a
   namesChile.txt file, which is the equivalent to the namesMap.txt file
   output by Xyce.  The program mapMerge.C is supposed to read these two files,
   and create a map between the chilespice and Xyce solution vectors.
 */

        if (matflag != 0 && ilast_time != 0)
        {
          sprintf(text, "x_vec%d", iterno);
          fp1      = fopen(text, "w");
          sprintf(text, "namesChile.txt");
          fp2      = fopen(text, "w");

          size     = SMPmatSize(ckt->CKTmatrix);
          icomplex = 0;

          for (i = 1; i <= size; i++)
            if (ckt->CKTirhs[i] != 0.0) icomplex = 1;

          for (i = 1; i <= size; i++)
          {
            if (icomplex == 1)
              fprintf(fp1, "%25.18e  %25.18e\n", ckt->CKTrhs[i],
                      ckt->CKTirhs[i]);
            else
              fprintf(fp1, "%25.18e\n", ckt->CKTrhs[i]);

            fprintf(fp2, "\t%d\t%s\n",i,CKTnodName(ckt,i) );
          }
          fclose(fp1);
          fclose(fp2);
        }

