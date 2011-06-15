#ifndef __DEBUGING_TOOLS_H__
#define __DEBUGING_TOOLS_H__

void print_sparse(MSKlidxt *ptrb, MSKlidxt *ptre, MSKidxt *sub, MSKrealt *val, int num_rows, int num_nonzeros) {
   int j;
   printf("ptrb:\n");
   for(j=0; j<num_rows; j++) {
      printf("%d ",ptrb[j]);
   }
   printf("\n");

   printf("ptre:\n");
   for(j=0; j<num_rows; j++) {
      printf("%d ",ptre[j]);
   }
   printf("\n");

   printf("sub:\n");
   for(j=0; j<num_nonzeros; j++) {
      printf("%d ",sub[j]);
   }
   printf("\n");

   printf("val:\n");
   for(j=0; j<num_nonzeros; j++) {
      printf("%f ",val[j]);
   }
   printf("\n");
}


void print_bounds(MSKboundkeye* bkx, MSKrealt* blx, MSKrealt* bux, MSKintt num_variables) {
   int j;
   printf("keys:\n");
   for(j=0; j<num_variables; j++) {
      printf("%d ",bkx[j]);
   }
   printf("\n");

   printf("lower bounds:\n");
   for(j=0; j<num_variables; j++) {
      printf("%f ",blx[j]);
   }
   printf("\n");

   printf("upper bounds:\n");
   for(j=0; j<num_variables; j++) {
      printf("%f ",bux[j]);
   }
   printf("\n");
}


void print_symmetric(MSKidxt* subi, MSKidxt* subj, MSKrealt *val, MSKintt num_nonzeros) {
   int j;
   printf("subi:\n");
   for(j=0; j<num_nonzeros; j++) {
      printf("%d ",subi[j]);
   }
   printf("\n");

   printf("subj:\n");
   for(j=0; j<num_nonzeros; j++) {
      printf("%d ",subj[j]);
   }
   printf("\n");

   printf("val:\n");
   for(j=0; j<num_nonzeros; j++) {
      printf("%f ",val[j]);
   }
   printf("\n");
}

#endif //  __DEBUGING_TOOLS_H__
