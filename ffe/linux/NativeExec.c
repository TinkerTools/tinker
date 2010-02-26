#include <stdlib.h>
#include <jni.h>

char *GetStringNativeChars(JNIEnv *env, jstring jstr);

JNIEXPORT jint JNICALL Java_ffe_core_FFEExec_NativeExec
   (JNIEnv *env, jobject obj, jstring jargv,
                              jstring jdir,
                              jstring jpath,
                              jstring jclasspath,
                              jstring jld_library_path) {
   char *argv;
   char *dir;
   char *path;
   char *classpath;
   char *ld_library_path;
   char *temp;
   char *current;
   int len;

   /* Turn Java strings into native strings */
   path = GetStringNativeChars(env, jpath);
   if (path == NULL) return -1;
   classpath = GetStringNativeChars(env, jclasspath);
   if (classpath == NULL) return -1;
   ld_library_path = GetStringNativeChars(env, jld_library_path);
   if (ld_library_path == NULL) return -1;
   dir = GetStringNativeChars(env, jdir);
   if (dir == NULL) return -1;
   argv = GetStringNativeChars(env, jargv);
   if (argv == NULL) return -1;

   /* Set the PATH */
   current = (char *) getenv("PATH");
   if (current != NULL) {
      temp = strstr(current, path);
      if (temp == NULL){
         len  = strlen(path) + strlen(current) + 2;
         temp = (char *) malloc(len);
         sprintf(temp, "%s:%s", path, current);
         setenv("PATH",temp,1);
         free(temp);
      }
   } else {
      setenv("PATH",path,1);
   }

   /* Set the CLASSPATH */
   current = (char *) getenv("CLASSPATH");
   if (current != NULL) {
      temp = strstr(current, classpath);
      if (temp == NULL) {
         len  = strlen(classpath) + strlen(current) + 2;
         temp = (char *) malloc(len);
         sprintf(temp, "%s:%s", classpath, current);
         setenv("CLASSPATH",temp,1);
         free(temp);
      }
   } else {
      setenv("CLASSPATH",classpath,1);
   }

   /* Set the LD_LIBRARY_PATH */
   current = (char *) getenv("LD_LIBRARY_PATH");
   if (current != NULL) {
      temp = strstr(current, ld_library_path);
      if (temp == NULL) {
         len  = strlen(ld_library_path) + strlen(current) + 2;
         temp = (char *) malloc(len);
         sprintf(temp, "%s:%s", ld_library_path, current);
         setenv("LD_LIBRARY_PATH",temp,1);
         free(temp);
      }
   } else {
      setenv("LD_LIBRARY_PATH",ld_library_path,1);
   }
   //printf("%s\n\n", ld_library_path);
   //printf("%s\n\n", path);
   //printf("%s\n\n", classpath);

   /* Change Directories and Execute the TINKER Command */
   chdir(dir);
   system(argv);

   free(path);
   free(classpath);
   free(ld_library_path);
   free(dir);
   free(argv);
   return 0;
}

char *GetStringNativeChars(JNIEnv *env, jstring jstr) {
   jbyteArray bytes = 0;
   jthrowable exc;
   jclass class;
   static jmethodID MID_String_getBytes = NULL;
   char *result = 0;

   if ((*env)->EnsureLocalCapacity(env, 2) < 0) {
      printf("Out of Memory\n");
      return 0; /* Out of Memory */
   }

   /* Do we have the cached Method ID Yet? */
   if (MID_String_getBytes == NULL) {
      class = (*env)->FindClass(env, "Ljava/lang/String;");
      if (class == 0) {
         printf("Could not find String Class\n");
         return 0; /* Could not Find the String Class */
      }
      MID_String_getBytes = (*env)->GetMethodID(env,
                                    class, "getBytes", "()[B");
      if (MID_String_getBytes == 0) {
         printf("Could not find getBytes\n");
         return 0;
      }
   }

   bytes = (*env)->CallObjectMethod(env, jstr, MID_String_getBytes);
   exc = (*env)->ExceptionOccurred(env);
   if (!exc) {
      jsize len = (*env)->GetArrayLength(env, bytes);
      result = (char *) malloc(len + 1);
      if (result == 0) {
         (*env)->DeleteLocalRef(env, bytes);
         return 0;
      }
      (*env)->GetByteArrayRegion(env,bytes,0,len,(jbyte *)result);
      result[len] = 0;
   } else {
      printf("Exception calling String->GetBytes\n");
      (*env)->DeleteLocalRef(env, exc);
   }
   (*env)->DeleteLocalRef(env, bytes);
   return result;
}
