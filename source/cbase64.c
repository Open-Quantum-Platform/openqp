#include <stddef.h>
#include <stdint.h>

static const char BASE64_TABLE[64] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
static const char BASE64_RTABLE[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 0, 0, 0, 63, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    0, 0, 0, 0, 0, 0, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

static const char BASE64_PAD = '=';

/*
 * @brief Encode `nbytes` bytes of `src` to base64 string at `dst`
 * @param src       input data
 * @param dst       output data
 * @param nbytes    number of bytes in `src`
 * @return number of bytes in encoded string
 */
int64_t
base64_encode(void *src, char *dst, long long int nbytes) {

    char *s = src;
    char *d = dst;

    long long int i;

    for (i = 0; i < nbytes - 2; i += 3 ) {
        *dst++ = BASE64_TABLE[                                       (unsigned char)s[i  ]>>2 ];
        *dst++ = BASE64_TABLE[ (unsigned char)(s[i  ] & 0x03) << 4 | (unsigned char)s[i+1]>>4 ];
        *dst++ = BASE64_TABLE[ (unsigned char)(s[i+1] & 0x0f) << 2 | (unsigned char)s[i+2]>>6 ];
        *dst++ = BASE64_TABLE[ (unsigned char) s[i+2] & 0x3f                                  ];
    }

/* remainder */
    switch (nbytes % 3) {
    case 0:
        break;
    case 1:
        *dst++ = BASE64_TABLE[ (unsigned char)                             s[i  ]>>2 ];
        *dst++ = BASE64_TABLE[ (unsigned char)(s[i  ] & 0b00000011) << 4             ];
        *dst++ = BASE64_PAD;
        *dst++ = BASE64_PAD;
        break;
    case 2:
        *dst++ = BASE64_TABLE[ (unsigned char)                             (unsigned char)s[i  ]>>2 ];
        *dst++ = BASE64_TABLE[ (unsigned char)(s[i  ] & 0b00000011) << 4 | (unsigned char)s[i+1]>>4 ];
        *dst++ = BASE64_TABLE[ (unsigned char)(s[i+1] & 0b00001111) << 2                            ];
        *dst++ = BASE64_PAD;
        break;
    }
    return dst - d;
}

/*
 * @brief Decode base64 string `src` to the byte array `dst`
 * @param src       base64-encoded input data
 * @param dst       output data
 * @return number of bytes readed from src, negative value means error on that position in input
 */
int64_t
base64_decode(char *src, char *dst) {

    int32_t ibuf;
    unsigned char itmp;
    char stmp;
    char obuf[3];

    int64_t nbytes = 0;

    int k = 0;
    int pad = 0;

    while (stmp = *src++) {
        nbytes++;
        switch ( stmp ) {
            case 'A'...'Z':
            case 'a'...'z':
            case '0'...'9':
            case '+':
            case '/':
                itmp = BASE64_RTABLE[stmp];
                break;
            case '=':
                itmp = 0;
                pad++;
                break;
            default:
                return -nbytes;
        }

        if (!pad) {
            ibuf = ibuf << 6 | (0x3f & itmp);
        } else {
            ibuf = ibuf << 6;
        }
        k = (k+1) % 4;

        if (!k) {
            ibuf = ibuf >> (8*pad);
            for (int j = 0; j < 3-pad; j++) {
                obuf[j] = ibuf & 0xff;
                ibuf = ibuf >> 8;
            }
            for (int j = 0; j < 3-pad; j++) {
                *dst++ = obuf[2-pad-j];
            }
            ibuf = 0;
        }
    }

/* remainder */
    if (k) {
        for (int j = k; j<4; j++) {
            ibuf = ibuf << 6;
        }

        for (int j = 0; j < 3; j++) {
            obuf[j] = ibuf & 0xff;
            ibuf = ibuf >> 8;
        }
        for (int j = 0; j < 3; j++) {
            *dst++ = obuf[2-j];
        }
    }

    return nbytes;

}
