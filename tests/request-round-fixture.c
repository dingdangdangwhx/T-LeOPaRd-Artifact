/* Public, synthetic input only. Not linked into a protocol binary. */
int fixture_client_cli(int argc, char **argv);
#define main fixture_client_cli
#include "../src/protocol/client.c"
#undef main
#include "fips202.h"

static uint64_t candidate;
static unsigned int draw;
void __wrap_randombytes(uint8_t *, size_t);

void __wrap_randombytes(uint8_t *out, size_t n) {
  uint8_t in[32] = "request-round-fixture-v1";
  if(n != 16 || draw >= 128) abort();
  /* Draw zero generates R and its opening; keep the statement fixed.
     Subsequent contexts use independent, candidate-specific test coins. */
  uint64_t domain = draw ? ((candidate << 8) | draw) : 0;
  for(size_t i=0;i<8;i++) in[24+i]=(uint8_t)(domain >> (8*i));
  shake128(out,n,in,sizeof(in));
  draw++;
}

int main(int argc, char **argv) {
  char *end;
  if(argc != 2 || !*argv[1]) return 1;
  candidate=strtoull(argv[1],&end,10);
  if(*end || candidate > UINT32_MAX) return 1;
  leopard_pp pp;
  leopard_pp_init(&pp);
  if(leopard_write_pp_file("fixture.pp",&pp)) return 1;
  char *args[]={"fixture","epoch-request","fixture.pp","6cb0da21981144bd9a457b80469d539e",
                "3","2","0,2","504f5052462d7461672d303030303100","5ac31980277104bed261",
                "fixture.req","fixture.st",NULL};
  int rc=fixture_client_cli(11,args);
  return rc || draw < 2;
}
