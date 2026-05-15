# GitHub workflows

The continuous integration tests (`tests.yml`) run with every push on a pull request.

The PyPI workflow (`publish_pypi.yml`) runs when a tag is pushed to GitHub, e.g.

```
git tag v1.2.3
git push origin v1.2.3
```

Make sure the git tag matches the STIsim version.