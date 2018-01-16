#include <decaf.hxx>
#include <decaf/eddsa.hxx>
#include <decaf/spongerng.hxx>
#include <stdio.h>

using namespace decaf;

void printbe(const SecureBuffer &b, bool clip_hi = false, const char *after="\n") {
    uint8_t mask = clip_hi ? 0x7F : 0xFF;
    printf("0x");
    for (ssize_t i=b.size()-1; i>=0; i--) {
        printf("%02x", b[i] & mask);
        mask = 0xFF;
    }
    printf("%s",after);
}

void describe(const IsoEd25519::Point &p) {
    printf("(");
    printbe(IsoEd25519::Scalar(4)*p,false,",");
    printbe(p.mul_by_half_cofactor_and_encode_like_eddsa(),true,")");
}

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;
    
    SpongeRng rng(Block("like_eddsa"),SpongeRng::DETERMINISTIC);
    
    printf("points = [\n");
    for (int i=0; i<100; i++) {
        if (i) printf(",\n");
        describe(IsoEd25519::Point(rng));
    }
    printf("\n]\n");
    
    printf("base = ");
    describe(IsoEd25519::Point::base());
    
    printf("\nbase_cof = ");
    printbe(IsoEd25519::Point::base().mul_by_half_cofactor_and_encode_like_eddsa());
    
    printf("\nbase_cof_448 = ");
    printbe(Ed448Goldilocks::Point::base().mul_by_half_cofactor_and_encode_like_eddsa());
    
    return 0;
}