while getopts 'q:o:' flag; do
  case "${flag}" in
    q)
      QUERY="${OPTARG}"
      ;;
    o)
      OUTPUT="${OPTARG}"
      ;;
    ?)
      echo "script usage: $(basename \$0) [-q query.fasta] [-o output_folder]" >&2
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

sed "s/ /_/g" $QUERY > $OUTPUT