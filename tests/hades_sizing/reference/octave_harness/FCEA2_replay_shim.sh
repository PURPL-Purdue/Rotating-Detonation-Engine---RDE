#!/bin/bash
# Replays committed CEA2 output only if the deck written by the MATLAB matches the committed deck.
read name
name=$(echo "$name" | tr -d '\r')
if cmp -s "$name.inp" "/tmp/hades_octave/committed/$name.inp"; then
  cp "/tmp/hades_octave/committed/$name.out" "$name.out"
  echo "FCEA2 shim: replayed committed $name.out" >&2
else
  echo "FCEA2 shim: deck $name.inp differs from committed deck; refusing" >&2
  diff "$name.inp" "/tmp/hades_octave/committed/$name.inp" >&2
  rm -f "$name.out"; exit 1
fi
