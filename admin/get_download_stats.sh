echo "GitHub release downloads"
curl -s -i https://api.github.com/repos/kusterlab/picked_group_fdr/releases -H "Accept: application/vnd.github.manifold-preview+json" | grep "\"name\"\|\"download_count\""
#https://pypistats.org/api/packages/simsi_transfer/recent

echo "PyPI downloads"
curl -s -i https://pypistats.org/api/packages/picked_group_fdr/recent -H "Accept: application/vnd.github.manifold-preview+json" | grep "\"last_day\""
