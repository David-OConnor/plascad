# This file sets up a Linux desktop entry, and moves the application to the home directory.

NAME_UPPER="PlasCAD"
NAME="plascad"

printf "Moving the ${NAME_UPPER} executable and icon to ~/${NAME}..."

chmod +x $NAME

if [ ! -d ~/${NAME} ]; then
  mkdir ~/${NAME}
fi

cp ${NAME} ~/${NAME}
cp icon.png ~/${NAME}/icon.png

# Update the desktop entry with the absolute path.
sed "s|~|$HOME|g" ${NAME}.desktop > ~/.local/share/applications/${NAME}.desktop

printf "\nComplete. You can launch ${NAME_UPPER} through the GUI, eg search \"${NAME_UPPER}\", and/or add to favorites.\n"