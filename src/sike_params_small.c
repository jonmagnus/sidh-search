#include <sike_params.h>

const sike_params_raw_t SIKEp33 = {
    .name = "SIKEp33",

    .p = "0x01037e1fff",

    .lA = "2",
    .eA = "13",
    .lB = "3",
    .eB = "12",
    
    .A=6,
    .B=1,

    .xQA0 = "0x7EDDCC08",
    .xQA1 = "0x2EECFC7F",
    .yQA0 = "0xCC97F758",
    .yQA1 = "0x803C8E58",
    .xPA0 = "0xC9DF8BA2",
    .xPA1 = "0xD7612D84",
    .yPA0 = "0x8C960D73",
    .yPA1 = "0x0E6F0B44",
    .xQB0 = "0x4C989623",
    .xQB1 = "0x044E0BC1",
    .yQB0 = "0x81148BAB",
    .yQB1 = "0x24899449",
    .xPB0 = "0x363492D9",
    .xPB1 = "0x5EF2B7A1",
    .yPB0 = "0x60B8E935",
    .yPB1 = "0x9EC02750",
    
    .crypto_bytes=4,
    .msg_bytes=4,
};
