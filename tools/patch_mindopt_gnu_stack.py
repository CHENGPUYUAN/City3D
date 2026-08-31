#!/usr/bin/env python3
"""Clear the PF_X bit of PT_GNU_STACK in MindOpt's shipped libraries.

MindOpt 2.3.0 ships libmindopt.so.2.3.0 with a RWE GNU_STACK flag, which
Ubuntu 26.04+ kernels refuse to load ("cannot enable executable stack").
Run this after every MindOpt (re)install:

    uv run tools/patch_mindopt_gnu_stack.py

It patches both copies (keeps a .bak next to each):
  ~/mindopt/<ver>/linux64-x86/lib/libmindopt.so.<ver>
  ~/mindopt/<ver>/linux64-x86/lib/python/mindoptpy/libmindopt.so.<ver>
"""
import glob
import os
import struct
import sys


def clear_px(path: str) -> bool:
    with open(path, "rb") as f:
        data = bytearray(f.read())
    e_phoff = struct.unpack_from("<Q", data, 0x20)[0]
    e_phentsize = struct.unpack_from("<H", data, 0x36)[0]
    e_phnum = struct.unpack_from("<H", data, 0x38)[0]
    patched = False
    for i in range(e_phnum):
        off = e_phoff + i * e_phentsize
        if struct.unpack_from("<I", data, off)[0] != 0x6474E551:  # PT_GNU_STACK
            continue
        flags = struct.unpack_from("<I", data, off + 4)[0]
        if flags & 1:
            bak = path + ".bak"
            if not os.path.exists(bak):
                with open(bak, "wb") as f:
                    f.write(data)
            struct.pack_into("<I", data, off + 4, flags & ~1)
            patched = True
            print(f"patched {path}: GNU_STACK {flags:#x} -> {flags & ~1:#x}")
    if patched:
        with open(path, "wb") as f:
            f.write(data)
    else:
        print(f"ok (already clean): {path}")
    return patched


def main() -> int:
    home = os.path.expanduser("~/mindopt")
    targets = glob.glob(f"{home}/*/linux64-x86/lib/libmindopt.so.*") + \
              glob.glob(f"{home}/*/linux64-x86/lib/python/mindoptpy/libmindopt.so.*")
    targets = [p for p in targets if not p.endswith(".bak")]
    if not targets:
        sys.exit(f"no libmindopt.so found under {home}; is MindOpt installed?")
    for t in targets:
        clear_px(t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
