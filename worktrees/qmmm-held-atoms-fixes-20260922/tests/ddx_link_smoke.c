#include <stdio.h>
#include <string.h>

#include "ddx.h"

int main(void) {
  char banner[2048];
  memset(banner, 0, sizeof(banner));
  ddx_get_banner(banner, (int)sizeof(banner));

  if (banner[0] == '\0') {
    fprintf(stderr, "ddX banner is empty\n");
    return 1;
  }

  printf("ddX link smoke test passed\n");
  return 0;
}
