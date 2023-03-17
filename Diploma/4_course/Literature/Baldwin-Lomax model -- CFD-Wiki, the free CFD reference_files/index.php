/** CSS placed here will be applied to all skins */
/* Infobox template style */

.infobox {
   border: 1px solid #aaa;
   background-color: #f9f9f9;
   color: black;
   margin-bottom: 0.5em;
   margin-left: 1em;
   padding: 0.2em;
   float: right;
   clear: right;
}
.infobox td,
.infobox th {
   vertical-align: top;
}
.infobox caption {
   font-size: larger;
   margin-left: inherit;
}
.infobox.bordered {
   border-collapse: collapse;
}
.infobox.bordered td,
.infobox.bordered th {
   border: 1px solid #aaa;
}
.infobox.bordered .borderless td,
.infobox.bordered .borderless th {
   border: 0;
}

.infobox.sisterproject {
   width: 20em;
   font-size: 90%;
}

@media print {
    .infobox.sisterproject {
        display: none;
    }
}

/* styles for bordered infobox with merged rows */
.infobox.bordered .mergedtoprow td,
.infobox.bordered .mergedtoprow th {
   border: 0;
   border-top: 1px solid #aaa;
   border-right: 1px solid #aaa;
}

.infobox.bordered .mergedrow td,
.infobox.bordered .mergedrow th {
   border: 0;
   border-right: 1px solid #aaa;
}