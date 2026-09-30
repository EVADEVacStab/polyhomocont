# polyhomocont

A Fortran library for finding all isolated complex solutions of square
systems of polynomial equations with homotopy continuation.

This package was written by Claude (Anthropic), including the source
code, the test programs and this README.

**Status:** under development. See `CLAUDE.md` for the roadmap.

## Building and testing

```sh
fpm build --profile release
fpm test --profile release
```

Always pass `--profile debug` or `--profile release` explicitly: without
it, fpm compiles with no flags at all, so none of the flags defined in
`fpm.toml` are applied.

## License

MIT, see [LICENSE](LICENSE).
