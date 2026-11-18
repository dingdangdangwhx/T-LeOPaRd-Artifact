#include "lazer.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

int main(void) {
  const unsigned int widths[]={12,23,24,29,30,31,32,34,40,52};
  INTVEC_T(input,11,2);
  INTVEC_T(decoded,11,2);
  INTVEC_T(prefix,1,2);
  INTVEC_T(prefix_dec,1,2);
  INT_T(mod,2);
  coder_state_t enc,dec;
  uint8_t bytes[1024];
  lazer_init();
  for(unsigned int w=0;w<sizeof(widths)/sizeof(widths[0]);w++) {
    int64_t a=INT64_C(1)<<(widths[w]+1);
    const int64_t values[]={-3*a-1,-a,-a/2-1,-a/2,-1,0,1,a/2-1,a/2,a,3*a+1};
    for(unsigned int j=0;j<11;j++) int_set_i64(intvec_get_elem(input,j),values[j]);
    for(unsigned int offset=0;offset<8;offset++) {
      memset(bytes,0,sizeof(bytes));
      coder_enc_begin(enc,bytes);
      if(offset) {
        int_set_i64(mod,INT64_C(1)<<offset);
        int_set_i64(intvec_get_elem(prefix,0),(INT64_C(1)<<offset)-1);
        coder_enc_urandom(enc,prefix,mod,offset);
      }
      coder_enc_grandom(enc,input,widths[w]); coder_enc_end(enc);
      coder_dec_begin(dec,bytes);
      if(offset) coder_dec_urandom(dec,prefix_dec,mod,offset);
      coder_dec_grandom(dec,decoded,widths[w]);
      if(!coder_dec_end(dec) || !intvec_eq(input,decoded) ||
         (offset && !intvec_eq(prefix,prefix_dec))) {
        fprintf(stderr,"Gaussian codec mismatch: width=%u offset=%u\n",widths[w],offset);
        return 1;
      }
    }
  }
  puts("Gaussian codec: widths 12..52, signs, boundaries, eight bit offsets PASS");
  return 0;
}
